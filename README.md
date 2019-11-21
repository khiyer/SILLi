# SILLi
A 1D Numerical Tool Quantifying the Thermal Effects of Sill Intrusions

### Introduction

SILLi is a user-friendly 1D FEM based tool which calculates the thermal effects of sill intrusions on the enclosing sedimentary stratigraphy and can be globally applied. Input data for the model is the present-day well log or sedimentary column with an Excel input file and includes rock parameters such as thermal conductivity, total organic carbon (TOC) content, porosity, and latent heats. The model accounts for sedimentation and burial based on a rate calculated by the sedimentary layer thickness and age. Erosion of the sedimentary column is also included to account for realistic basin evolution. Multiple sills can be emplaced within the system with varying ages.

The model output includes the thermal evolution of the sedimentary column through time, and the changes that take place following sill emplacement such as TOC changes, thermal maturity, and the amount of organic and carbonate-derived CO<sub>2</sub>. The TOC and vitrinite results can be readily benchmarked within the tool to present-day values measured in the sedimentary column. This allows the user to determine the conditions required to obtain results that match observables and leads to a better understanding of metamorphic processes in sedimentary basins.

### Governing Equation

The governing equation is the 1D heat diffusion equation:

![Governing Equation](img/image001.png)
 
Where <I>c<sub>peff</sub></I> is the effective rock heat capacity that includes the latent heat of crystallization for the sill and <I>H</I> accounts for the latent heat of dehydration in the surrounding sedimentary rock.

### Examples
- sills in the Karoo basin (v1.0) [Karoo Basin](/tutorials/karoo.md)
- sills in the Vøring basin (v1.0) [Vøring Basin](/tutorials/utgard.md)
- igneous pluton without basin evoultion (v1.0) [igneous intrusion into a pluton](/tutorials/pluton.md) 
- sills in a coal bed methane play in Botswana (v1.1) (Bulguroglu and Milkov, 2020)

### Further Reading
A full description of the model with examples can be found [here](https://www.geosci-model-dev.net/11/43/2018/).
