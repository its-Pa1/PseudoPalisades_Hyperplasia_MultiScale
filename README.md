# **Repository Overview – PseudoPalisades_Hyperplasia_MultiScale**  

Welcome to this GitHub repository, which houses the code associated with the article **"A Flux-Limited Model for Glioma Patterning with Hypoxia-Induced Angiogenesis"**, available at [MDPI](https://www.mdpi.com/2073-8994/12/11/1870). The results are also included in **Chapter 3** of the thesis, accessible at [KLUEDO](https://kluedo.ub.rptu.de/frontdoor/index/index/docId/6573).  

This repository primarily focuses on the **thesis results**, which provide a **more detailed and comprehensive analysis** compared to the article. The thesis covers everything from the paper while incorporating additional insights and extended discussions.  

## **Key Features**  

This repository provides implementations for:  

- **Full model simulations** of glioma patterning  
- **Comparison between different models** of tumor growth and diffusion  
- **Pattern formation analysis** in the tumor microenvironment  

## **Project Structure**  

### **Scripts**  
The main script, **`main_all.m`**, runs all the individual scripts mentioned below. The results are presented in **Chapter 3** of the thesis (*Diss_Kumar_Pawan.pdf* in the parent directory).  

### **Full_model_I_II** (Figures 3.2 & 3.3)  
Contains codes and results for **experiments 1 and 2** described in **Chapter 3** of the thesis.  

### **Difference_minimal** (Figure 3.4)  
Contains codes and results comparing the **full model** with the **old model** (minimal model), specifically looking at glioma density differences within **experiment 2**.  

### **Difference_self_diffusion** (Figure 3.5)  
Includes codes and simulations comparing the **full model** with a version **without limited flux** for glioma self-diffusion (*γ₂ = 0, denominator of ph_taxis = 1*) within **experiment 2**.  

### **Pattern** (Figure 3.6 – First & Second Row)  
This directory contains **1D codes and videos** for the **full system in dimensional form**, along with **pattern formation results**.  

### **Pattern_difference** (Figure 3.7)  
Includes codes and results showing the **pattern differences** between the **old and new models** within **experiment 3 (1D case)**.  

### **Pattern_diff_W_L_flux** (Figure 3.6 – Third to Sixth Row)  
Contains simulations of pattern formation for:  
1. The model **without self-diffusion**  
2. The **difference** between the full model and the model **without self-diffusion**  

---

Each subdirectory contains a separate **README** file explaining its contents in detail.  

Feel free to explore the directories to dive deeper into specific aspects of the project. If you have any questions or need further clarification, don't hesitate to reach out.  

## **Author**  
**Pawan Kumar**: [@its-Pa1](https://github.com/its-Pa1)  

© 2025, Pawan Kumar. All Rights Reserved.  
