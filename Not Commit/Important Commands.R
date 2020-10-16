#Creating package license
usethis::use_mit_license("Rodolfo Pelinson")

#This is to link the package to a github repository
usethis::use_git()
usethis::use_github()

#Creating a README file
usethis::use_readme_rmd()

#Creating functions
#The name of the function should alwais end with "_function"
usethis::use_r("<NAME>_function")

#Creating data to the package
#add the data objects inside the parentesys
usethis::use_data()
#You should also create a "data.R" folder inside the "R" folder. You will add your data information there.


#Typing this, devtools will generate appropriate documentation
devtools::document()



usethis::use_data(com_SS2_SS3_abundance,
                  isolation_SS2_SS3,
                  fish_SS2_SS3,
                  SS_SS2_SS3,
                  ID_SS2_SS3,
                  All)


usethis::use_data(Trait_SS2_SS3)
usethis::use_data(com_SS2_SS3)

com_SS2_SS3_non_predators_abundance <- rowSums(com_SS2_SS3_non_predators)
com_SS2_SS3_predators_abundance <- rowSums(com_SS2_SS3_predators)

usethis::use_data(com_SS2_SS3_non_predators_abundance,com_SS2_SS3_predators_abundance)

usethis::use_data()
fish_isolation_SS2_SS3 <- fish_isolation

