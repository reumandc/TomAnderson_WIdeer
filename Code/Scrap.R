<!--checkpoint package for CRAN-available packages
```{r checkpoint_chunk, echo=F, warning=F, message=F, results="hide"}
library(checkpoint)
if (!dir.exists("./.checkpoint/")){
  dir.create("./.checkpoint/")
}
checkpoint("2018-04-20",checkpointLocation = "./")
```
-->
  

Things to do:
D - checkpoint and install_github need to play together, somehow
T - set.seed in appropriate places make it exactly reproducible
D - put pedagog fig code in
DT - use cacheing, which means sorting out dependencies
sort out warnings, don't just supress them
ready to add data to repo?