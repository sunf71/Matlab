**************************************************************
**************************************************************
Single Path Video Stabilization
Version 1.0
Written by Shuaicheng Liu (liuyangmao@gmail.com), date: 2014.Nov.15

Implementation of the basic single path 2D video stabilization
1. estimate single homographies between neighboring frames
2. smooth a single camera path

This code is directly translated from our C++ version. 
This code does not contain adaptive path smoothing, thus cannot handle camera motions such as quick rotation and zooming.
As the purpose is to demonstrate the basic path smoothing strategy, we hope it can lower the barrier for the implementation of the 2D video stabilization. 
It only estimates a single path, and adopts the basic path smoothing strategy. However, it has the potential to stabilize many challenge videos captured by consumer-grade devices. 




If you use/adapt our code in your work, you need to appropriately cite our SIGGRAPH 2013 paper.

This code is for academic purpose only. Not for commercial/industrial activities.





**************************************************************
**************************************************************
Usage: 

run demo.m in Matlab


Acknowledgement: Thanks Heng Guo for testing and video I/O. 



