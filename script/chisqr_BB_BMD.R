BB_users = c(215, 248)
BB_nonuser = c(660, 618)
table <- rbind(BB_users, BB_nonuser)
name(table)<- c('BMD_low', 'BMD_normal')
Chi <- chisq.test(table)
Chi
