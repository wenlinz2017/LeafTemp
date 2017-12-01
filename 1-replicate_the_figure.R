vp3 <- function(Tmpt) #a function to estimate VPD unite: kpa source: Campbell&Norman Env. biophysics
{return(0.611*exp((17.502*Tmpt)/(Tmpt+240.97)))}

# a function to estimate the slope of the saturation vapor pressure function
Delta.fun <- function(x)
{return(17.502*240.47*vp3(x)/(240.97+x)^2)}


leaf.T <- function(R_ni,g_vs) # a function to estimate leaf temperature
{
T_a <- 30
u <- 1
RH <- 0.2

gamma <- 6.66E-4 #thermodynamic psychrometer constant
p_a <- 101.3
c_p <- 29.3 # specific heat of air at constant pressure, J/mol/K

g_r <- 7E-6*(T_a+273.15)^2-0.0021*(T_a+273.15)+0.2085;g_r # radiative conductance

e_s <- vp3(T_a);e_s # saturated vapor pressure 

Delta <- Delta.fun(T_a);Delta # slope of vapor pressure function
s <- Delta/p_a # slope of saturation mole fraction function
D <- e_s*(1-RH);D


d <- 1:100/1000 # leaf dimension

g_Ha <- 0.135*sqrt(u/d);g_Ha # example 14.1 used g_Ha <- 1.4*0.135*sqrt(u/d)
g_va <- 0.147*sqrt(u/d);g_va # example 14.1 used g_va <- 1.4*0.147*sqrt(u/d)
g_Hr <- g_Ha+g_r;g_Hr

g_v <- g_vs*g_va/(g_vs+g_va);g_v
gamma_star <- gamma *g_Hr/g_v;gamma_star

T_l <- T_a + gamma_star/(gamma_star+s)*(R_ni/(g_Hr*c_p)-D/p_a/gamma_star);T_l
return(T_l-T_a)
}

dat1 <- leaf.T(300,0.4)
dat2 <- leaf.T(300,0.01)
dat3 <- leaf.T(0,0.4)
dat4 <- leaf.T(0,0.01)

par(mar=c(4.5,4.5,1,1))
plot(log10(1:100),dat1,type="l",col="mistyrose3",ylim=c(-10,20),xaxt="n",
lwd=2,ylab=expression(T[l]~"-"~T[a]~"("~degree~C~")"),xlab="Leaf dimensions (mm)")
axis(1,at=log10(c(1:10,seq(20,100,by=10))),labels=F)
axis(1,at=c(0,1,2),labels=c(1,10,100))
lines(log10(1:100),dat2,col="indianred3")
lines(log10(1:100),dat3,col="violet")
lines(log10(1:100),dat4,col="lightpink3")
abline(h=0,col=8,lty=3)
legend("topleft",c("Rni=300,gvs=0.4","Rni=300,gvs=0.01","Rni=0,gvs=0.4",
      "Rni=0,gvs=0.01"),lty=1,col=c("mistyrose3","indianred3",
      "violet","lightpink3"))

