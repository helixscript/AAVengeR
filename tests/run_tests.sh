sudo docker run -it --mount type=bind,source=/data,dst=/data/ aavenger_v135_env bash 

source /root/.bashrc

source activate bioinformatic_dependencies

aavenger_test_dir="/data/AAVengeR/tests"

Rscript $aavenger_test_dir/check_config_file.R -h

# Rscript $aavenger_test_dir/check_config_file.R --run_type build --overwrite

Rscript $aavenger_test_dir/check_config_file.R --run_type test

# Rscript $aavenger_test_dir/core_module_integration_test.R --run_type build --output_dir /data --number_of_cpus 4

Rscript $aavenger_test_dir/core_module_integration_test.R --run_type test --output_dir /data/fake_test2 --number_of_cpus 32