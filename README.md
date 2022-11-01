# independencepvalue-project

This R package is created using literate programming with the  [litr](https://github.com/jacobbien/litr-project/tree/main/litr) R package.  Please see [independencepvalue](independencepvalue) for the generated R package itself.

## Code for generating the `independencepvalue` package

After cloning this repo, please use

```r
remotes::install_github("jacobbien/litr-project", subdir = "litr")
litr::render("create-independencepvalue.Rmd")
```

This will create [create-independencepvalue.html](https://htmlpreview.github.io/?https://github.com/ArkajyotiSaha/independencepvalue-project/blob/main/create-independencepvalue.html) and the package directory [independencepvalue](independencepvalue).
