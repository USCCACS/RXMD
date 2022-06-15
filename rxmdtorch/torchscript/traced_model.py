import torch
import torchvision
from torch import nn

# An instance of your model.
model = torchvision.models.resnet18()

# An example input you would normally provide to your model's forward() method.
example = torch.rand(1, 3, 224, 224)

# Use torch.jit.trace to generate a torch.jit.ScriptModule via tracing.
traced_script_module = torch.jit.trace(model, example)
traced_script_module.save("traced_resnet_model.pt")

# another example with branch
class MyModule(torch.nn.Module):
    def __init__(self, N, M):
        super(MyModule, self).__init__()
        self.weight = torch.nn.Parameter(torch.rand(N, M))

    def forward(self, input):
        if input.sum() > 0:
          output = self.weight.mv(input)
        else:
          output = self.weight + input
        return output

class MyNeuralNetwork(nn.Module):
    def __init__(self, input_size, output_size):
        super(MyNeuralNetwork, self).__init__()
        self.nn = nn.Sequential(
                nn.Linear(input_size, 20),
                nn.ReLU(),
                nn.Linear(20, 20), 
                nn.ReLU(),
                nn.Linear(20, output_size)
        )

    def forward(self, x):
        return self.nn(x)

my_module = MyModule(10,20)
sm = torch.jit.script(my_module)
sm.save("my_module_model.pt")

feature_size=21*3
my_nn = MyNeuralNetwork(feature_size,1)
example = torch.rand(feature_size)
traced_script_module = torch.jit.trace(my_nn, example)
traced_script_module.save("traced_nn_model.pt")

print(my_nn(example))
