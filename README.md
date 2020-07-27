# Improving Solar Cell Designs using Convolutional Neural Networks

**Sumit Bhattacharya , Devanshu Arya , ........................,Rajat Mani Thomas, Deepak Kumar Gupta**

#### Abstract
Topology Optimization has been a very popular method for optimizing physical structures. This has also been extended to the domains of Heat Flow and Solar Cell design problems. But, the performance is dependent on proper selection of parameters to be optimized. Better designs with better performance have been obtained with the use of Deep Learning in structural optimization. Metallization designs for Solar Cells obtained by topology optimization tend to give good efficiency. However,very limited use of deep learning has been made in optimizing Solar Cell metallization patterns. In this paper, we present the use of Convolutional Neural Networks without any training data to obtain more robust and efficient designs.


#### Method
Solar Cell Metallization problem is a non-linear optimization task. It can be formulated as a minimization task ,with density variable lying between 0 and 1, and no volume constraint. Those densities are then used to compute the loss ( in our case, the loss is negative of output power), which required to be minimized. Conventionally, the density variables are optimized using algorithms like MMA and OC. 

We have tried to reparameterize the domain of variables to be optimized. We have employed a Convolutional Neural Network (CNN) , to output the density variables. Here, we try to optimize the weights and the biases of the CNN , to produce such values of density variables that will give us maximum efficiency (minimum loss).

This is a schematic diagram of our optimization process

   **---------------------------------------------------------------------------------------> Forward Pass**
   
   ![Test Image1](https://github.com/BhattacharyaSumit/deeptop_sol/blob/master/Figs/Screenshot%20from%202020-07-26%2012-04-00.png)     ![test 2](https://github.com/BhattacharyaSumit/deeptop_sol/blob/master/Figs/Screenshot%20from%202020-07-26%2012-34-36.png)
    **<---------------------------------------------------------------------------------------------------Gradients**

