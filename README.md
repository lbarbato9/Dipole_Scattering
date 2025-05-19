# Electromagnetic Field radiated by a Dipole in the presence of a layered spherical structure: model of a system head-high permittivity material in Magnetic Resonance

Ultra-High Field Magnetic Resonance Imaging (UHF MRI) refers to MRI systems
operating at magnetic field strengths of 7 Tesla and above. In particular, UHF-MRI
provides SNR values higher than those obtained in typical 3T and 1.5T MRI.
Despite this MRI-UFH also brings along some challenging issues that need to be
addressed. The main disadvantages are the field inhomogeneity due to comparable
wavelengths with tissue sizes and an increase in SAR. For general solving of such kind
of issues, very expensive and computationally complex hardware is used, which limits
the large-scale development.
High-permittivity materials represent a great alternative; in fact, the obtained re-
sults have been excellent.To predict the permittivity value that maximizes the field
within the subject being examined, a new analytical electromagnetic model of scattering in spherical geometry is developed, which allows replacing the traditional numerical
model.
The new model provides the evaluation of the field inside a multi-layer sphere,
representing a fast tool that allow to make prediction before considering experiment on
the human tissue directly.
Generally, the field in MRI is irradiated by coils, which can provide very good results
at low field intensities. It is even possible to increase the SNR simply by increasing the
number of loops. However, when the intensities reach 7T or more, this cannot be done,
since in that case the increase in loops does not influence the SNR.
Electric dipoles, on the other hand, can be an interesting solution to explore,especially
with high frequencies. That is why the analytical model has been implemented so far
considering the field irradiated by electric dipoles

1. Conversion from Mie Scattering to Inward-Outward Representation
The Mie scattering has been studied extensively and an explanation of the mathematics behind this problem can be found in the following book : "Julius Adams Stratton. Electromagnetic theory, volume 33. John Wiley & Sons,
2007". Starting from this theory the problem was reformulated considering the field as a superposition of progressive and regressive waves : "Giuseppe Ruello and Riccardo Lattanzi. Scattering from spheres: A new look into
an old problem. Electronics, 10(2):216, 2021".

2. Our contribution fits in the middle. Because the analytical method allows us to model the phenomenon of fields inside spherical geometries leaving out the Mie coefficients and using engineering parameters such as impedances and reflection coefficients.Considering this conversion we introduced an electric dipole as a field source and implemented everything on MATLAB. We used as startint point the following papers : "Jorge R Zurita-S´anchez. Anapole due to a tangentially polarized dipole near a sub-
wavelength sphere: Inhibition and enhancement of spontaneous emission. Physical
Review A, 104(5):053524, 2021." & "R Ruppin. Decay of an excited molecule near a small metal sphere. The Journal
of Chemical Physics, 76(4):1681–1684, 1982."

3. The convalidation of the results has been done comparing our simulation with the ones coming from a numerical solver.



  












