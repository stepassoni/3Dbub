# 3Dbub
3Dbub is a comprehensive repository that provides AI-based methodology and tools for image segmentation of two-phase flows in corrugated channels. This repository includes trained weights for the segmentation models, annotated images used for training and MATLAB scripts for bubble volume reconstruction and void fraction estimation.

## Features
- Pre-trained weights for image segmentation models based on U-Net architecture tailored to two-phase flow in corrugated channels.
- MATLAB scripts for accurate volume reconstruction and void fraction estimation.
- Annotated images used for training.

## Contents
- `matlab_scripts/`: Collection of MATLAB scripts for volume reconstruction and void fraction estimation.
- `model_weights/`: Directory containing the pre-trained weights for the U-net segmentation model ("3Dbub.caffemodel.h5", once extracted) and the model architecture ("3Dbub.modeldef.h5").
- `training_images/`: The images used for model training with annotations embedded (need to be opened in Fiji to see the annotations).

## Prerequisites
This repository requires the use of the U-Net plugin [[1]](#1) implemented in Fiji [[2]](#2). A complete guide on installation and usage is available [HERE](https://lmb.informatik.uni-freiburg.de/resources/opensource/unet/).

## Usage
To segment images using the pre-trained weights, you need to load the model file "3Dbub.modeldef.h5" and its weights "3Dbub.caffemodel.h5" in the U-Net plugin available in Fiji. Model has been finetuned to segment images of two-phase flow in chevron-type corrugated channels after background subtraction as detailed in the methodology described [HERE](https://doi.org/10.1016/j.ijmultiphaseflow.2024.104871). The input files are a set of frames from high-speed videos in the form of multi-page tiff files. The output from the segmentation plugin should be converted to 8-bit gray and then can be read into the Matlab scripts for further processing.

## Citation
Please cite the following work if you use this project's software, data, or methodologies in your own research:

>Stefano Passoni, Riccardo Mereu, Matteo Bucci,
>Integrating machine learning and image processing for void fraction estimation in two-phase flow through corrugated channels,
>_International Journal of Multiphase Flow_,
>Volume 177,
>2024,
>104871,
>ISSN 0301-9322,
>https://doi.org/10.1016/j.ijmultiphaseflow.2024.104871.

```bibtex
@article{PASSONI2024,
title = {Integrating machine learning and image processing for void fraction estimation in two-phase flow through corrugated channels},
journal = {International Journal of Multiphase Flow},
volume = {177},
pages = {104871},
year = {2024},
issn = {0301-9322},
doi = {https://doi.org/10.1016/j.ijmultiphaseflow.2024.104871},
url = {https://www.sciencedirect.com/science/article/pii/S0301932224001484},
author = {Stefano Passoni and Riccardo Mereu and Matteo Bucci},
}
```
## License
This project is licensed under the Creative Commons Attribution 4.0 license. See the [LICENSE](./LICENSE.txt) file for details.

## References
<a id="1">[1]</a> 
Falk, T., Mai, D., Bensch, R. et al. U-Net: deep learning for cell counting, detection, and morphometry. Nat Methods 16, 67–70 (2019). [doi](https://doi.org/10.1038/s41592-018-0261-2)

<a id="2">[2]</a> 
Schindelin, J., Arganda-Carreras, I., Frise, E., Kaynig, V., Longair, M., Pietzsch, T., … Cardona, A. (2012). Fiji: an open-source platform for biological-image analysis. Nature Methods, 9(7), 676–682. [doi](https://doi.org/10.1038/nmeth.2019)
