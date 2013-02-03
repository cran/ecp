 getWithin <- function (alpha_, X_) 
.Call("getWithin", alpha_, X_, PACKAGE = "ecp")

 getBetween <- function (alpha_, X_, Y_) 
.Call("getBetween", alpha_, X_, Y_, PACKAGE = "ecp")

 splitPointC <- function (s_, e_, D_, min_size_) 
.Call("splitPointC", s_, e_, D_, min_size_, PACKAGE = "ecp")

