# usethis::create_github_token()  
# gitcreds::gitcreds_set()       
# usethis::use_git_config(user.name = "Zhentao20", user.email = "zyu18@ad.unc.edu")

getwd()
list.files("/Users/zhentaoyu/Desktop/Project_1", all.files = TRUE)


# usethis::use_git()
usethis::proj_set("/Users/zhentaoyu/Desktop/Project_1")
usethis::proj_sitrep()
usethis::use_git_ignore(c(".Rproj.user", ".Rhistory", ".RData", ".Ruserdata"))

usethis::use_github()

install.packages("gert")  
gert::git_add(".")
gert::git_commit("Initial commit")

usethis::use_github(protocol = "https")