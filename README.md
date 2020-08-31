# TI
trajectory inference methods

测试、分解、打包一篇优秀的NBT的trajectory inference方法

他们总结了很多TI的方法，有很多测评，以及直接可运行的代码，非常实用！！！

[A comparison of single-cell trajectory inference methods](https://www.nature.com/articles/s41587-019-0071-9)

[Benchmarking trajectory inference methods](https://github.com/dynverse/dynbenchmark/)

[dynverse](https://dynverse.org/)

代码仓库，可直接改用：[dynverse: benchmarking, constructing and interpreting single-cell trajectories](https://github.com/dynverse)
- definition.yml，所有的参数以及默认值
- run.R，原始的R wrapper
- Dockerfile，docker运行代码

存在的问题：
- 他们wrapper写得有些问题，必须装了docker或Singularity才能运行
- 集群上没有docker或Singularity
- 个人Mac的16G内存不够docker或Singularity用
- 无法探究wrapper里面的机制，做个性化修改

解决方案：
- 好在他们提供了wrapper的所有代码
- 我把他们的wrapper改成我能利用的格式就好
- 同时还能熟悉算法的细节
