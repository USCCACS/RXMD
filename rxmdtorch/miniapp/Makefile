# Polaris
polaris:
	CC -O2 main.cpp -L../build -L/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-Wl,-R/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch -lstdc++

# Sunspot
sunspot:
	icpx main.cpp -O2 \
	-L../build \
	-L/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages/torch/lib \
	-Wl,-R/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages/torch/lib \
	-L/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages/intel_extension_for_pytorch/lib \
	-Wl,-R/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages/intel_extension_for_pytorch/lib \
	-ltorch_cpu -ltorch  -lc10 -lrxmdtorch -lintel-ext-pt-gpu -lintel-ext-pt-cpu


clean:
	rm a.out
