## 0.1. **tsBNgen: A Python Library to Generate Time Series Data Based on an Arbitrary Bayesian Network Structure**

[Description](#Description)

[Citation](Citaton)

[Features](#Features)

[Instruction](#Instruction)

[License](#License)

[Support](#Support)

----

### **Description**

### tsBNgen is a Python package to generate time series data based on an arbitrary Bayesian Network Structures. 
---
### **Citation**

 ### If you find this package useful or if you use it in your research or work please consider citing it as follows:
----
### **Features**

 - It handles discrete nodes, continous nodes and hybrid (Mixture of discrete and continuous) network.

 - It uses multinomila distribution for the discrete nodes and Gaussian distribution for the continuous nodes.

 - It handles arbitrary Bayesian network structure.

 - It supports loopback 1 (node at time t-1 can be connected to any node at time t) and loopback 2.

 - The code can be modified easily to handle arbitrary long loopbacks.
---

### **Instruction**

 To run this code either clone this repo or use the package distribution in PyPI using the following commands:

```python
pip install tsBNgen
```

 Then Run through the set of examples in 
 
 > **Time_Series_Generation_Examples.ipynb**

For more information on how to use the package please visit the followings:

1. Original paper 
2. Documentation in PDF available in this repo.

### **License**

This is repo is under the MIT liecense.

### **Support**

If you have any question or comment, either open a new issue for this repo or please email me at manitadayon@ucla.edu.













