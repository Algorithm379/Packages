#Package function set for testing
#Call data stored with package or follow vignette instructions for your own test analysis

#----------------------------------------------------------------

---
title: "Fiscore Package: Effective protein structural data visualisation and exploration "
author: "Austė Kanapeckaitė"
date: "8/10/2021"
output: html_document

---




##  Preprocessing PDB files----------------------------------------------------------------


pdb_path<- system.file("extdata", "6kz5.pdb", package="Fiscore")


pdb_df<-PDB_process(pdb_path)
pdb_df


##  Preparing PDB files----------------------------------------------------------------


pdb_path<-system.file("extdata", "6kz5_A.pdb", package="Fiscore")
pdb_df<-PDB_prepare(pdb_path)
head(pdb_df)


## Dihedral angle distribution plot (Ramachandran plot with densities)



phi_psi_plot(pdb_df)



## Dihedral angle distribution plot



phi_psi_bar_plot(pdb_df)




##  Dihedral angle distribution and secondary structure interactive plot



phi_psi_interactive(pdb_df)



##  Interactive scaled B-factor distribution


B_plot_normalised(pdb_df)


## Interactive 3D plot visualising dihedral angles and scaled B-factor distribution


phi_psi_3D(pdb_df)



## Interactive Fi-score plot visualising score distribution----------------------------------------------------------------


Fi_score_plot(pdb_df)

##   Calculate Fi-score for an individual region


Fi_score_region(pdb_df,50,70)


##  Fi-score plot visualising secondary structure elements

Fiscore_secondary(pdb_df)


## Hydrophobicity plot exploration----------------------------------------------------------------


hydrophobicity_plot(pdb_df,window = 9,weight = 25,model = "linear")


hydrophobicity_plot(pdb_df,window = 9,weight = 25,model = "exponential")



## Gaussian Mixture Models for structural feature prediction and classification----------------------------------------------------------------


df<-cluster_ID(pdb_df)

#predefined option implementation
df<-cluster_ID(pdb_df,clusters = 5, modelNames = "VVI")

## Project summary through density plots----------------------------------------------------------------
  

density_plots(pdb_df)

#including GMM output
density_plots(pdb_df, df)


