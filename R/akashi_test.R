akashi_test <- function(contingency_file=NULL){
    contingency_data <- read.table(contingency_file,sep=',')
    colnames(contingency_data) <- c('a','b','c','d')
    a <- contingency_data$a
    b <- contingency_data$b
    c <- contingency_data$c
    d <- contingency_data$d
    n <- a + b + c + d
    E_a <- (a + b)*(a + c)/n
    V_a <- (a + b)*(a + c)*(b + d)*(c + d)/(n^2*(n-1))
    Z <- (sum(a) - sum(E_a))/sqrt(sum(V_a))
    p <- 1 - pnorm(Z)
    odd_ratio <- sum(a*d/n)/sum(b*c/n)
    con_op=sum(a)
    con_nop=sum(c)
    var_op=sum(b)
    var_nop=sum(d)
    results <- list(Z=Z, 
                    p=p, 
                    odd_ratio=odd_ratio, 
                    conserved_optimal_sites=con_op,
                    conserved_non_optimal_sites=con_nop,
                    variable_optimal_sites=var_op,
                    variable_non_optimal_sites=var_nop,
                    con_var_ratio=(con_op+con_nop)/(var_op+var_nop)
                    )
    return(results)
}

