all:
	make polaris

# Polaris
polaris:
	CC -O2 main.cpp -L../build -L/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-Wl,-R/soft/datascience/conda/2022-07-19/pytorch/torch/lib \
	-ltorch_cpu -ltorch_cuda -lc10 -lrxmdtorch -lstdc++

# Sunspot
#sitepkg=/soft/datascience/aurora_models_frameworks-2023.0/lib/python3.9/site-packages
sitepkg=/home/knomura/ipex
sunspot:
	icpx main.cpp -O2 \
	-L../build \
	-L$(sitepkg)/torch/lib -Wl,-R$(sitepkg)/torch/lib \
	-L$(sitepkg)/intel_extension_for_pytorch/lib -Wl,-R$(sitepkg)/intel_extension_for_pytorch/lib \
	-ltorch_cpu -ltorch  -lc10 -lrxmdtorch -lintel-ext-pt-gpu -lintel-ext-pt-cpu

clean:
	rm a.out
