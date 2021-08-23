#call data stored in a package for testing---------------------------------------------------------------

pdb_path<- system.file("extdata", "6kz5.pdb", package="Fiscore")

#split files---------------------------------------------------------------
pdb_df<-PDB_process(pdb_path)

#output
pdb_df


#explore PDB input---------------------------------------------------------------
pdb_path<-system.file("extdata", "6kz5_A.pdb", package="Fiscore")
pdb_df<-PDB_prepare(pdb_path)

head(pdb_df)




#Plot---------------------------------------------------------------
phi_psi_plot(pdb_df)
phi_psi_bar_plot(pdb_df)
phi_psi_interactive(pdb_df)
phi_psi_3D(pdb_df)
Fi_score_plot(pdb_df)
Fi_score_region(pdb_df,50,70)
Fiscore_secondary(pdb_df)

#Hydrophobicity exploration---------------------------------------------------------------
hydrophobicity_plot(pdb_df,window = 9,weight = 25,model = "linear")
hydrophobicity_plot(pdb_df,window = 9,weight = 25,model = "exponential")

#Machine learning GMM implementation---------------------------------------------------------------

df<-cluster_ID(pdb_df)

#select specific model (not automatic processing implementation)
df<-cluster_ID(pdb_df,clusters = 5, modelNames = "VVI")


#Density plots and data summary---------------------------------------------------------------

density_plots(pdb_df)

#Machine learning combination with density plots
density_plots(pdb_df, df)

