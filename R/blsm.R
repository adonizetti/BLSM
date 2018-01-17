proc.crr=function(Z,Z0){
  #' @title Procrustean corresponding positions
  #' @description Given a set of starting coordinates, the function returns the Procrustean Transform of the initial points that minimizes 
  #' the sum of squared positional difference from a set of reference coordinates. The (Euclidean) distances between a candidate 
  #' configuration and the reference are evaluated by considering the couples of corresponding points. 
  #' 
  #' @param Z set of initial coordinates to be transformed
  #' @param Z0 set of reference coordinates
  #' 
  #' @return Set of coordinates minimizing the distance between the initial configuration and the reference one
  #' 
  #' @export
  
  Z=t(t(Z)-colMeans(Z)+colMeans(Z0))
  
  A=t(Z)%*%(Z0%*%t(Z0))%*%Z
  eA=eigen(A,symmetric=T)
  Ahalf=eA$vec%*%diag(sqrt(eA$val))%*%t(eA$vec)
  
  t(t(Z0)%*%Z%*%solve(Ahalf)%*%t(Z))  
}

gg_colors = function(n) {
  #' @title Color palette definition
  #' @description  Color palette definition
  #' @param n number of colors
  #' 
  #' @return Color palette
  
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

estimate_latent_positions = function (Y,W,
                                      procrustean = TRUE, 
                                      k=2,
                                      alpha=2,
                                      nscan=8*10^5, burn_in=5*10^5, odens=10^3,
                                      zdelta=10, z_norm_prior_mu=0, z_norm_prior_sd=10,
                                      adelta=.3, a_exp_prior_a=1, a_exp_prior_b=1,
                                      dynamic_plot = FALSE, dynamic_circles = FALSE,
                                      ...){
  #' @title Core BLSM simulation
  #' @description Run a simulation to obtain the positions of the network nodes in the latent space for each sampled iteration
  #' 
  #' @param Y Adjacency matrix of the network
  #' @param W Weight matrix of the network
  #' @param k Space dimensionality
  #' @param procrustean Boolean to include/exclude (TRUE/FALSE) the Procrustean Transform step in the algorithm. Set TRUE by default.
  #' @param alpha Starting value of the \code{alpha} parameter
  #' @param nscan Number of iterations
  #' @param burn_in Burn-in value (starting iterations to be discared)
  #' @param odens Thinning: only 1 iteration every \code{odens} will be sampled and stored in the output
  #' @param zdelta Standard deviation of the Gaussian proposal for latent positions
  #' @param z_norm_prior_mu Mean of the Gaussian prior distribution for latent positions 
  #' @param z_norm_prior_sd Standard deviation of the Gaussian prior distribution for latent positions
  #' @param adelta The uniform proposal for \code{alpha} is defined on the \eqn{[-adelta,+adelta]} interval
  #' @param a_exp_prior_a Shape parameter of the Gamma prior distribution for \code{alpha}. As the value is usually set to 1 the prior is an exponential distribution.
  #' @param a_exp_prior_b Rate parameter of the Gamma prior distribution for \code{alpha}. 
  #' @param dynamic_plot Boolean to plot dinamically the simulated positions (one update every \code{odens} iterations)
  #' @param dynamic_circles Boolean to add circles of radius \code{alpha} to the dynamic plots
  #' @param {\dots} Additional parameters that can be passed to \link[BLSM]{plot_latent_positions}
  #' 
  #' @return Returns a "BLSM object" (\code{blsm_obj}), i.e. a list containing:
  #' \itemize{
  #' \item \code{Alpha} {\code{alpha} values from the sampled iterations}
  #' \item \code{Likelihood}{Log-likelihood values from the sampled iterations}
  #' \item \code{Iterations}{Latent space coordinates from the sampled iterations. Latent positions are stored in a
  #' 3D array whose dimensions are given by (1: number of nodes, 2: space dimensionality, 3: number of iterations).
  #' In the non-Procrustean framework the latent distances are given instead of the positions: another 3D array is returned, whose dimensions
  #' are given by (1: number of nodes, 2: number of nodes, 3: number of iterations). The command needed in order to get the average values over the iterations for
  #' either the positions or the distances is \code{rowMeans(blsm_obj$Iterations, dims=2)}.}
  #' \item \code{StartingPositions}{Latent space coordinates right after the initialization step. In the non-Procrustean framework starting distances are given instead.}
  #' \item \code{Matrix}{Original matrices of the network (adjacency and weights)}
  #' \item \code{Parameters}{List of parameters specified during the call to \link[BLSM]{estimate_latent_positions}}
  #' }
  #' 
  #' @export
  
  if (k==3 & dynamic_plot){
    require(rgl)
  }
  params = list(procrustean=procrustean, 
                k=k,
                alpha=alpha,
                nscan=nscan, burn_in=burn_in, odens=odens,
                zdelta=zdelta, z_norm_prior_mu=z_norm_prior_mu, z_norm_prior_sd=z_norm_prior_sd,
                adelta=adelta, a_exp_prior_a=a_exp_prior_a, a_exp_prior_b=a_exp_prior_b)
  
  nscan = nscan + burn_in
  it_cont = 1
  escape_flag=FALSE
  create_window_flag=dynamic_plot
  
  rem = which(rowMeans(Y)==0)
  if(length(rem)>0) {
    Y=Y[-rem,-rem]
    W=W[-rem,-rem]
  }
  n=dim(Y)[1]   
  my_colors=gg_colors(n)
  cc=(Y>0)+0
  D=dst(cc)
  Z=cmdscale(D, k)
  Z[,1]=Z[,1]-mean(Z[,1])
  Z[,2]=Z[,2]-mean(Z[,2])
  
  tmp_opt=c(alpha,c(Z))
  tmp_opt=optim(tmp_opt,mlpY,Y=Y,W=W, method="SANN")$par 
  tmp_opt=optim(tmp_opt,mlpY,Y=Y,W=W,method="Nelder-Mead")$par
  
  tmp_opt=tmp_opt*2/(tmp_opt[1])  
  alpha=tmp_opt[1]
  
  Z_Proc = matrix(tmp_opt[-1],nrow=n,ncol=k)
  Z = Z_Proc
  lpz = lpz_dist(Z)
  
  acc_a=0
  acc_z=0 
  Alpha=alpha       
  Lik=lpY(Y,lpz,alpha, W) 
  
  inputs = list(Adjacency = Y, Weight = W)
  
  if (procrustean){
    blsm_obj=list(Alpha=rep(NA, (nscan-burn_in)/odens), 
                  Likelihood=rep(NA, (nscan-burn_in)/odens),
                  Iterations=array(NA,dim=c(n,k,(nscan-burn_in)/odens)), 
                  StartingPositions=Z,
                  Matrix=inputs, 
                  Parameters=params)
    
    avg_Z_est = Z
    for(ns in 1:nscan){
      tmp = Z_up(Y,Z,W,alpha,zdelta, z_norm_prior_mu, z_norm_prior_sd)
      if(any(tmp!=Z)){
        acc_z=acc_z+sum(tmp!=Z)/(2*n*odens)
        Z[,1] = Z[,1]-mean(Z[,1])
        Z[,2] = Z[,2]-mean(Z[,2])
        tryCatch({
          Z = proc.crr(tmp,Z_Proc)
        }, error = function(e) {
          escape_flag<<-TRUE
          message("The matrix used to compute the Procrustean transformation is singular. \nIf changing the parameters doesn't solve the issue, please try to lower the space dimensionality.")
          graphics.off()
        }
        )
      }
      if (escape_flag){
        return(blsm_obj)
      }
      lpz = lpz_dist(Z)
      tmp = alpha_up(Y,lpz,W,alpha,adelta,a_exp_prior_a,a_exp_prior_b)
      if(tmp!=alpha) {
        acc_a=acc_a+1/odens
        alpha=tmp
      }
      if (ns%%odens==0){
        
        if(ns>burn_in){
          if (create_window_flag){
            if (k==2) {
              x11(xpos=0,ypos=0)
            } else if (k==3) {
              par3d(windowRect = c(50, 50, 800, 800))
            } else {
              message("Error: plot cannot be displayed since space dimensionality is bigger than 3.")
              dynamic_plot=FALSE
            }
            create_window_flag = FALSE
          }
          blsm_obj$Alpha[it_cont]= alpha
          lik = lpY(Y, lpz, alpha, W)
          blsm_obj$Likelihood[it_cont] = lik
          cat(ns-burn_in,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          blsm_obj$Iterations[,,it_cont] = Z
          avg_Z_est= ((it_cont-1)*avg_Z_est+Z)/it_cont 
          it_cont = it_cont+1 
          
          if (dynamic_plot){
            plot_latent_positions(blsm_obj, circles_2D = dynamic_circles, ...)
          }
        } else {
          if (ns==odens){
            cat("\nBeginning burn-in period...\n\n")
          }
          lik = lpY(Y, lpz, alpha, W)
          cat(ns,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          if (ns==burn_in){
            cat("\nBurn-in period ended.\n\nBeginning simulation...\n\n")
          }
        }
      }
    }
  } else {
    blsm_obj=list(Alpha=rep(NA, (nscan-burn_in)/odens), 
                  Likelihood=rep(NA, (nscan-burn_in)/odens),
                  Iterations=array(NA,dim=c(n,n,(nscan-burn_in)/odens)),
                  StartingDistances=-lpz_dist(Z),
                  Matrix=inputs, 
                  Parameters=params)
    
    for(ns in 1:nscan){
      Z_tmp=Z_up(Y,Z,W,alpha,zdelta, z_norm_prior_mu, z_norm_prior_sd)  
      if(any(Z_tmp!=Z)){
        acc_z=acc_z+sum(Z_tmp!=Z)/(2*n*odens)
        Z=Z_tmp
      }
      lpz = lpz_dist(Z)
      tmp=alpha_up(Y,lpz,W,alpha,adelta,a_exp_prior_a,a_exp_prior_b)
      if(tmp!=alpha) {
        acc_a=acc_a+1/odens
        alpha=tmp
      }
      if (ns%%odens==0){
        if( ns > burn_in){
          blsm_obj$Alpha[it_cont]= alpha
          lik = lpY(Y, lpz, alpha, W)
          blsm_obj$Likelihood[it_cont] = lik
          cat(ns-burn_in,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          blsm_obj$Iterations[,,(ns-burn_in)/odens] = -lpz
        } else {
          if (ns==odens){
            cat("\nBeginning burn-in period...\n\n")
          }
          lik = lpY(Y, lpz, alpha, W)
          cat(ns,acc_a,acc_z,alpha,lik,"\n")
          acc_z=acc_a=0
          if (ns==burn_in){
            cat("\nBurn-in period ended.\n\nBeginning simulation...\n\n")
          }
        }
      }
    }
    
  }
  cat("\nThe simulation ended successfully.")
  return(blsm_obj)
}

plot_traceplots_acf = function (blsm_obj, chosen_node=1,  coordinate=1, chosen_pair=c(1,2)){
  #' @title BLSM traceplots and ACF
  #' @description Traceplots and autocorrelation functions for the \code{alpha} parameter and a selected node (or pair of nodes in the non-Procrustean framework)
  #' 
  #' @param blsm_obj Blsm object obtained through \link[BLSM]{estimate_latent_positions}
  #' @param chosen_node Specified node for traceplot and autocorrelation function (Procrustean framework)
  #' @param coordinate Specified coordinate dimension from the n-dimensional latent space
  #' @param chosen_pair Specified pair of nodes for traceplot and autocorrelation function (non-Procrustean framework)
  #' 
  #' @export
  
  par(mfrow=c(2,2))
  if (blsm_obj$Parameters$procrustean==TRUE){
    plot(blsm_obj$Alpha,type="l",main="Alpha")
    plot(blsm_obj$Iterations[chosen_node,coordinate,],type="l", ylab=paste0("Coordinate ",coordinate), main = paste0("Node ",chosen_node))
    plot(acf(blsm_obj$Alpha,plot=F),main="Alpha")
    plot(acf(blsm_obj$Iterations[chosen_node,coordinate,],plot=F),main= paste0("Node ",chosen_node))
  } else {
    plot(blsm_obj$Alpha,type="l",main="Alpha")
    plot(blsm_obj$Iterations[chosen_pair[1],chosen_pair[2],],type="l", ylab="Euclidean distance", main = paste0("Nodes ",chosen_pair[1], " and ", chosen_pair[2]))
    plot(acf(blsm_obj$Alpha,plot=F),main="Alpha")
    plot(acf(blsm_obj$Iterations[chosen_pair[1],chosen_pair[2],],plot=F),main= paste0("Nodes ",chosen_pair[1], " and ", chosen_pair[2]))
  }
}

plot_latent_positions = function(blsm_obj, colors, points_size=0.1, labels_point_size=5, labels_point_color="yellow", 
                                 labels_text_size=1, labels_text_color="blue", circles_2D = FALSE){
  #' @title Base BLSM plot function
  #' @description Plot latent positions from a Procrustean simulation
  #' 
  #' @param blsm_obj Blsm object obtained through \link[BLSM]{estimate_latent_positions}
  #' @param colors Colors of the simulated coordinate points in the latent space
  #' @param points_size Size of the coordinate points
  #' @param labels_point_size Size of the label points
  #' @param labels_point_color Color of the label points
  #' @param labels_text_size Text size in the label points
  #' @param labels_text_color Text color in the label points
  #' @param circles_2D Plot circles of radius \code{alpha} (see the model's main variables) centered around the label points
  #' 
  #' @export
  
  n=dim(blsm_obj$Iterations)[1]
  if (missing(colors)) {
    colors=gg_colors(n)
  }
  
  if (blsm_obj$Parameters$procrustean==TRUE){
    if (dim(blsm_obj$Iterations)[2]==2){
      par(mfrow=c(1,1))
      par(mar=c(1,1,1,1))
      par(mgp=c(2,1,0))
      dev.hold()
      plot(blsm_obj$Iterations[,1,],blsm_obj$Iterations[,2,],pch=20,cex=points_size, col=colors, xlab="", ylab="", xaxt="n", yaxt="n")
      
      avg_Z_est=rowMeans(blsm_obj$Iterations, dims=2, na.rm = TRUE)
      
      for(i in 1:(n-1)){
        for(j in (i+1):n){
          lines(avg_Z_est[c(i,j),1], avg_Z_est[c(i,j),2] ,lty=blsm_obj$Matrix$Adjacency[i,j], col="blue", lwd =2)
        }
      }
      
      points(avg_Z_est[,1],avg_Z_est[,2],xaxt="n",yaxt="n",xlab="",ylab="", col=labels_point_color,pch=20,cex=labels_point_size)
      text(avg_Z_est[,1],avg_Z_est[,2],labels(blsm_obj$Matrix$Adjacency)[[1]],col=labels_text_color,cex=labels_text_size)
      
      if (circles_2D){symbols(avg_Z_est, circles = rep(mean(Alpha),n), add=TRUE, fg = colors, inches=F)}
      dev.flush()
    } else if (dim(blsm_obj$Iterations)[2]==3){
      par(mfrow=c(1,1))
      par(mar=c(3,3,1,1))
      par(mgp=c(2,1,0))
      
      plot3d(blsm_obj$Iterations[,1,],blsm_obj$Iterations[,2,],blsm_obj$Iterations[,3,], size=points_size*10, 
             col=colors, xlab = "", ylab="", zlab="")
      rgl.viewpoint(theta=30, phi=10, fov=30)
      
      if (missing(avg_Z_est)){
        avg_Z_est=rowMeans(blsm_obj$Iterations, dims=2, na.rm = TRUE)
      }
      for(i in 1:(n-1)){ 
        for(j in (i+1):n){
          if (blsm_obj$Matrix$Adjacency[i,j]){
            lines3d( avg_Z_est[c(i,j),], col="blue", lwd =2)
          }
        }
      }
      
      points3d(avg_Z_est[,1],avg_Z_est[,2],avg_Z_est[,3],xaxt="n",yaxt="n",xlab="",ylab="", col=labels_point_color,size=labels_point_size*10)
      
      text3d(avg_Z_est[,1],avg_Z_est[,2],avg_Z_est[,3],labels(blsm_obj$Matrix$Adjacency)[[1]],col=labels_text_color,cex=labels_text_size)
    } else {
      message("Error: plot cannot be displayed since space dimensionality is bigger than 3.")
    }
  } else {
    message("Sampled latent positions are not available in the non-Procrustean framework.\n")
  }
}
