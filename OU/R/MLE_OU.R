
 require (stats4)
 require (sde)




 dcOU <- function (x , t , x0 , theta , log = FALSE ){
   Ex <- theta[1] / theta [2]+( x0 - theta[1] / theta[2]) * exp ( - theta[2] * t)
   Vx <- theta[3]^2 * (1 - exp (-2 * theta[2] *t )) / (2 * theta[2])
   dnorm (x , mean = Ex , sd = sqrt (Vx), log = log )
 }
 OU.lik <- function ( theta1 , theta2 , theta3 ){
   n <- length (X)
   dt <- deltat (X)- sum (dcOU(X [2: n ], dt , X[1:( n -1)] , c(theta1 , theta2 , theta3 ) ,log = TRUE ))
 }



 set.seed (123)
 X <- sde.sim ( model = "OU" , theta = c (3 ,1 ,2) , N =1000 , delta =1)
 mle(OU.lik , start = list(theta1 =1 , theta2 =0.5 , theta3 =1) ,method = "L-BFGS-B" , lower = c( - Inf ,0 ,0)) -> fit
 summary ( fit )

  mle ( minuslogl = OU . lik , start = list ( theta1 = 1 , theta2 = 0.5 ,theta3 = 1) , method = "L - BFGS -B " , lower = c (- Inf , 0 , 0)



dcOU <- function (x , t , x0 , theta , log = FALSE ){
  Ex <- theta [1] / theta [2]+( x0 - theta [1] / theta [2]) * exp ( - theta [2] * t)
  Vx <- theta [3]^2 * (1 - exp ( -2 * theta [2] *t )) / (2 * theta [2])
  dnorm (x , mean = Ex , sd = sqrt ( Vx ), log = log )
}
OU . lik <- function ( theta1 , theta2 , theta3 ){
  n <- length ( X)
  dt <- deltat ( X)
  - sum ( dcOU (X [2: n ], dt , X [1:( n -1)] , c ( theta1 , theta2 , theta3 ) , log = TRUE ))
}
