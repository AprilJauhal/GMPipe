#EXAMPLE RUN: snakemake --directory test_dir --jobs 999 --wait-for-files -w 500 --cluster 'bsub -J snakemake -n 1 -o Pipe_snakemake.stdout -e Pipe_snakemake.stderr'
# Make sure to run GMPipe_start.sh first 

test_ORF, = glob_wildcards("storage/hmm/pass/{test_ORF}.fa")

rule ALL:
	input:
		"out/PASSING_SEQUENCES.fa",
		"storage/opt_ML/best_model.txt"
		
rule MUSCLE_PROFILES:
	input:
		"storage/hmm/pass/{each_ORF}.fa"
	output:
		"storage/ML/{each_ORF}/{each_ORF}.afa"
	shell:
		"source GMPipeline_userinput.ctl; "
		"${{MUSCLE}} -profile -in1 storage/reference.afa -in2 {input} -out {output}"

rule IQTREE_unconstr:
	input:
		"storage/ML/{each_ORF}/{each_ORF}.afa"
	output:
		treefile="storage/ML/{each_ORF}/{each_ORF}_unconstr{rep,\d+}.treefile",
		iq_file="storage/ML/{each_ORF}/{each_ORF}_unconstr{rep,\d+}.iqtree"
	shell:
		"source GMPipeline_userinput.ctl; "
		"BEST_MODEL=$(cat storage/opt_ML/best_model.txt); "
		"${{IQTREE}} -s {input} -m ${{BEST_MODEL}} -redo -nt AUTO -ntmax 1 -allnni -pre $(dirname {output.iq_file})/$(basename {output.iq_file} .iqtree)"

		
rule IQTREE_constr:
	input:
		"storage/ML/{each_ORF}/{each_ORF}.afa"
	output:
		treefile="storage/ML/{each_ORF}/{each_ORF}_constr{rep,\d+}.treefile",
		iq_file="storage/ML/{each_ORF}/{each_ORF}_constr{rep,\d+}.iqtree"
	shell:
		"source GMPipeline_userinput.ctl; "
		"BEST_MODEL=$(cat storage/opt_ML/best_model.txt); "
		"${{IQTREE}} -s {input} -m ${{BEST_MODEL}} -redo -g storage/opt_ML/main.constr -nt AUTO -ntmax 1 -allnni -pre $(dirname {output.iq_file})/$(basename {output.iq_file} .iqtree)"

rule CONCAT:
	input:
		"storage/ML/{each_ORF}/{each_ORF}_unconstr01.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr02.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr03.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr04.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr05.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr06.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr07.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr08.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr09.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_unconstr10.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr01.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr02.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr03.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr04.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr05.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr06.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr07.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr08.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr09.treefile",
		"storage/ML/{each_ORF}/{each_ORF}_constr10.treefile"
        
	output:
		"storage/ML/concat/{each_ORF}/{each_ORF}_comparison.treelst"
	shell:
		"cat {input} > {output} "
		
rule IQTREE_COMPARE:
	input:
		"storage/ML/concat/{each_ORF}/{each_ORF}_comparison.treelst"
	output:
		"storage/ML/concat/{each_ORF}/{each_ORF}_comparison.iqtree"
	shell:
		"source GMPipeline_userinput.ctl; "
		"BEST_MODEL=$(cat storage/opt_ML/best_model.txt); "
		"ORF_comparison_name=$(basename {output} .iqtree); "
		"ORF_name=${{ORF_comparison_name%_comparison}}; "
		"${{IQTREE}} -s storage/ML/${{ORF_name}}/${{ORF_name}}.afa -m ${{BEST_MODEL}} -redo -z {input} -n 0 -zb 10000 -au -pre $(dirname {output})/${{ORF_comparison_name}}"
		
rule FINALIZATION:
	input:  
		expand("storage/ML/concat/{each_ORF}/{each_ORF}_comparison.iqtree", each_ORF=test_ORF)
	output:
		"out/PASSING_SEQUENCES.fa"
	shell:
		"source GMPipeline_userinput.ctl; "
		"bash ${{SCRIPT_PATH}}/GMPipe_results.bat GMPipeline_userinput.ctl "
