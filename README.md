# Prediction-of-knee-joint-kinematics

This thesis focuses on the creation of a mathematical model of the knee joint capable of predicting the internal 
kinematics of the joint and the load on soft tissues during motion. A mathematical model was constructed based 
on the joint’s geometry, mechanical properties of individual structures, and calculation of external loads. The input
to the model is the angle of flexion in the knee joint, and the output is motion around three axes (three rotations 
and three translations) and the loading of joint structures. The resulting movements are consistent with experimentally 
measured knee joint kinematics during the stance phase of walking. The results of structural loading are in qualitative 
agreement with the literature.


Instruction:\
In input.py set model of ligaments (linear or non-linear), mechanical properties of knee issues and optimization 
parametrs.\
Then run main.py.\
All results  will be save in folder 'Results_linear' or 'Results_non-linear'.\
For video output run video.py.\
To compare the result with linear and non-linear ligaments run results.py.
