# EasyBrainSurfacePlots 0.0.1

This repository contains some modified functions based on the packages `brainspace`, `neuromaps` to help people who are unfamiliar with existing neuroscience tools in generating nice surface brain plots using only a few lines of code, starting from a vector of real values.

![Parcellated surface with colormaps](https://github.com/andresantoro/EasyBrainSurfacePlots/blob/main/easybrainsurfaceplot-multiview.jpg)


# Example code
With just a few lines of code, you can generate a picture similar to the one shown above.
To do so, clone this repository and run the Jupyter notebook  `01_example_surface_plots.ipynb`.
The notebook requires the following Python packages: `brainspace`, `neuromaps`, and `svgutils`. You can install them directly using `pip`, for example, `pip install brainspace`.

```
import sys
sys.path.append('.')
from utils_neuromaps_brain import *
import numpy as np
import matplotlib.pyplot as plt

node_strength=np.random.randint(0,100,size=100)
N=100
fig=normal_view(node_strength,edges=True,exp_form=False)
```


The code has been tested with the following versions of the python packages:
- brainspace==0.1.10 
- neuromaps==0.0.2
- matplotlib==3.6.3
- svgutils==0.3.2


## License

This project is licensed under the GNU GPLv3 license. Please see the [LICENSE](LICENSE) file for more details.

