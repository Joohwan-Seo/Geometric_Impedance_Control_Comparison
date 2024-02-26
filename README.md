# Geometric-Impedance-Control
Tested in Matlab R2021a and R2021b

Written by Joohwan Seo, Ph.D. Student in Mechanical Engineering, UC Berkeley.

## Authors' comment
These are quite raw, unorganized files. I hope everyone can get some of the insights, and please forgive me for the nasty codes.

## Fixes & Update

## Main files
`main_geo_discrete_bullo.m` is a GIC-1 Controller, a Frobenius-norm based, or Lie-group based control design.
`main_geo_discrete_log_map.m` is a GIC-2 Controller, a Lie-algebra based control design.
`plotter_geo_log_comp.m` for visualizing the comparison between GIC-1 and GIC-2 controller.\
`potential_comparison.m` for visualizing the comparison result for SO(2) case, i.e., Figure 2.\
### Note
Comparison results on `plotter2.m` are not presented in the paper since the two control laws are different, thus unable to do a fair comparison.

`plotter_geo_comp.m` is for comparison between intuitive geometric impedance control and GIC-v1. When the gains become different, they start to behave differently.\
When gains are scalar in translational and rotational ways, both controllers are identical.

## Codes directly implemented from:
https://github.com/Joohwan-Seo/Geometric-Impedance-Control-Public

## Objectives in Main files
Change `obj` variable to `tracking`, `tracking2`, `regulation`, and `regulation2`.\
`tracking`: sinusoidal trajectory tracking \
`tracking2`: 3rd-order smooth polynomial trajectory tracking both in translational and rotational (presented in paper)\
`regulation`: regulation task (presented in paper)\
`regulation2`: Multi-point regulation, can be considered as step-input case.

## In "trajectory_generator" folder
`desired_trajectory.m` gives the trajectory utilized in `tracking` objective.\
`desired_trajectory2.m` gives the trajectory utilized in `tracking2` objective.\
`desired_trajectory_regulation.m` gives the desired setpoint utilized in `regulation` objective.\
`desired_trajectory_regulation2.m` gives the desired setpoints utilized in `regulation2` objective.\
`trajectory_calculator.m` is for the implementation of 3rd-order polynomial smooth trajectory generation.\

## "function_generator" folder
`RTB matlab` (Robotics Toobox Matlab) is needed to run the code. \
Build the UR5e robot model and associated dynamic parameter matrices as well as Jacobian matrices.

## "sub_direct" folder
Dynamic parameter matrices built from function_generator are saved here. Some miscellaneous functions are also defined here.

## "results" folder
The simulation result data generated by the main files are saved here.

## Accepted and will be presented at
American Control Conference (ACC) 2024, Toronto, Canada \
``A Comparison Between Lie Group- and Lie Algebra- Based Potential Functions for Geometric Impedance Control''

Arxiv submitted version:
https://arxiv.org/abs/2401.13190

## Bibtex Citation
@article{seo2024comparison,\
  title={A Comparison Between Lie Group-and Lie Algebra-Based Potential Functions for Geometric Impedance Control},\
  author={Seo, Joohwan and Prakash, Nikhil Potu Surya and Choi, Jongeun and Horowitz, Roberto},\
  journal={arXiv preprint arXiv:2401.13190},\
  year={2024}\
}
