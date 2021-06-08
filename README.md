In this project, we provide the Matlab code of the point cloud coarse registration, which is implemented via 2D line features. 
The main advantage of this algorithm is its high computation efficiency. Everyone is welcome to use the code for research work, but not for commerce.
If you use the code, please cite my paper
(Wuyong Tao, Xianghong Hua, Zhiping Chen and Pengju Tian. (2020). Fast and automatic registration of terrestrial point clouds using 2D line features. Remote Sensing, 12(8): 1283-1298.)


In this project, three files are provided. The "line2D" file is used to extract 2D line feature from point clouds. 
The "partialoverlapICP3" file is used to perform the ICP (iterative closest point) algorithm, which is a popular fine registration method, after our algorithm. 
The "point cloud registration based on 2D line features" file is used to perform the coarse registration algorithm in our paper.

Before you carry out our algorithm, you need to calculate the point cloud resolution (pr). 
