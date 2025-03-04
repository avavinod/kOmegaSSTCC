# kOmegaSSTCC
This is an attempt to add curvature correction (CC) to the k-Omega SST turbulence model based on the paper by Smirnov and Mentner "Sensitization of the SST Turbulence Model to Rotation and Curvature by Applying the Spalart-Shur Correction term".

Motivation: It is well knownt that cyclone separator simulations employing the normal k-Omega SST and other such two equation models fail to predict the correct velocity distributions because of strong swirl and streamline curvature. The curvature correction term is a way of enabling the two-equation models to more accurately capture such flow behaviour without resorting to RSTM models.

How does it work: Uses empirical relations to modify the turbulent kinetic energy production term. 

How is the model implemented: I created a copy of the kOmegaSSTBase model under turbulenceModels/Base and added it to the RAS models folder under kOmegaSSTCC folder. The more elegant way would have been to modify the kOmegaSSTBase model with additional CC terms that can be switched on and off depending on the model employed. But in order to not break anything too much, the base model was copied.

IS IT WORKING: Not yet. Currently, I am trying to debug run-time some errors that are plaguing the implementation. Once it is working, I will update this page.