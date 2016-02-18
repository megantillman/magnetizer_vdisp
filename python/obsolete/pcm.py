""" this code compares output from multiple simulation runs """
import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MultipleLocator, FormatStrFormatter,MaxNLocator,NullFormatter,FixedLocator, AutoMinorLocator
from matplotlib import rc
from rcm import * #read in simulation output


pi=np.pi
sqrt=np.sqrt
atan=np.arctan

nt=420
ntfrac=0.8
print 'nt=%d' %nt
t1_Gyr=0.8
t2_Gyr=1.2
t1=t1_Gyr/t0_Gyr
t2=t2_Gyr/t0_Gyr
print t1_Gyr, t2_Gyr

plotdirec=s1+'plots/'

os.chdir(plotdirec)

r_kpc=         np.zeros((nruns,ngmax,     cols))

for i in range(nruns):
  for j in range(int(ngals[i])):
    r_kpc[i,j,:]=r[i,j,:]*r_disk_kpc[i,j]

ts_H=          ts_h/ts_l
ts_H[np.isnan(ts_H)]=0.
ts_Phi=        ts_Om*ts_tau
ts_V=          ts_Uz/ts_v
ts_V[np.isnan(ts_V)]=0.
ts_q=         -ts_G/ts_Om
ts_pkin_rad=  -sqrt(2./pi/ts_q)/ts_H    
ts_pkin_rad[np.isnan(ts_pkin_rad)]=0.
ts_pkin_deg=   180./pi*ts_pkin_rad
ts_Etot=       np.ndarray.sum(ts_Br**2+ts_Bp**2+ts_Bzmod**2,axis=3)  #total magnetic energy
ts_Btot=       sqrt(ts_Etot)


ts_t_Gyr=      np.zeros((nruns,ngmax,rows     ))
ts_Br_mkG=     np.zeros((nruns,ngmax,rows,cols))
ts_Bp_mkG=     np.zeros((nruns,ngmax,rows,cols))
ts_Bzmod_mkG=  np.zeros((nruns,ngmax,rows,cols))
ts_Btot_mkG=   np.zeros((nruns,ngmax,rows     ))
ts_Uphi_kms=   np.zeros((nruns,ngmax,rows,cols))
ts_alpk_kms=   np.zeros((nruns,ngmax,rows,cols))
ts_alpm_kms=   np.zeros((nruns,ngmax,rows,cols))
ts_alp=        np.zeros((nruns,ngmax,rows,cols))
ts_alp_kms =   np.zeros((nruns,ngmax,rows,cols))
ts_h_kpc=      np.zeros((nruns,ngmax,rows,cols))
ts_l_kpc=      np.zeros((nruns,ngmax,rows,cols))
ts_v_kms=      np.zeros((nruns,ngmax,rows,cols))
ts_n_cm3=      np.zeros((nruns,ngmax,rows,cols))
ts_Uz_kms=     np.zeros((nruns,ngmax,rows,cols))
ts_Om_Gyr=     np.zeros((nruns,ngmax,rows,cols))
ts_G_Gyr=      np.zeros((nruns,ngmax,rows,cols))
ts_etat_cm2s=  np.zeros((nruns,ngmax,rows,cols))
ts_tau_Gyr=    np.zeros((nruns,ngmax,rows,cols))
ts_Beq_mkG=    np.zeros((nruns,ngmax,rows,cols))
ts_Uphi_kms=   np.zeros((nruns,ngmax,rows,cols))

for i in range(nruns):
  for j in range(int(ngals[i])):
    ts_t_Gyr    [i,j,:]=  ts_t     [i,j,:] *t0_Gyr [i,j]
    ts_Br_mkG   [i,j,:]=  ts_Br    [i,j,:] *B0_mkG [i,j]
    ts_Bp_mkG   [i,j,:]=  ts_Bp    [i,j,:] *B0_mkG [i,j]
    ts_Bzmod_mkG[i,j,:]=  ts_Bzmod [i,j,:] *B0_mkG [i,j]
    ts_Btot_mkG [i,j,:]=  ts_Btot  [i,j,:] *B0_mkG [i,j]
    ts_alpk_kms [i,j,:]=  ts_alpk  [i,j,:] *h0_kpc [i,j]/t0_kpcskm[i,j]
    ts_alpm_kms [i,j,:]=  ts_alpm  [i,j,:] *h0_kpc [i,j]/t0_kpcskm[i,j]
    ts_alp      [i,j,:]=  ts_alpk  [i,j,:] +ts_alpm[i,j]
    ts_alp_kms  [i,j,:]=  ts_alp   [i,j,:] *h0_kpc [i,j]/t0_kpcskm[i,j]
    ts_h_kpc    [i,j,:]=  ts_h     [i,j,:] *h0_kpc [i,j]
    ts_l_kpc    [i,j,:]=  ts_l     [i,j,:] *h0_kpc [i,j]
    ts_v_kms    [i,j,:]=  ts_v     [i,j,:] *h0_kpc [i,j]/t0_kpcskm[i,j]
    ts_n_cm3    [i,j,:]=  ts_n     [i,j,:] *n0_cm3 [i,j]
    ts_Uz_kms   [i,j,:]=  ts_Uz    [i,j,:] *h0_kpc [i,j]/t0_kpcskm[i,j]
    ts_Om_Gyr   [i,j,:]=  ts_Om    [i,j,:] /t0_Gyr [i,j]
    ts_G_Gyr    [i,j,:]=  ts_G     [i,j,:] /t0_Gyr [i,j]
    ts_etat_cm2s[i,j,:]=  1./3*ts_l[i,j,:]*ts_v[i,j,:]*etat0_cm2s[i,j]/1.e26
    ts_tau_Gyr  [i,j,:]=  ts_l[i,j,:]/ts_v[i,j,:]*t0_Gyr[i,j]
    ts_Beq_mkG  [i,j,:]=  sqrt(4*pi*ts_n[i,j,:])*ts_v[i,j,:]*B0_mkG[i,j]
    for k in range(cols):
      ts_Uphi_kms [i,j,:,k]=ts_Om[i,j,:,k]*r[i,j,k]*r_disk_kpc[i,j]/t0_kpcskm[i,j]
ts_tau_Myr=     ts_tau_Gyr*1000
ts_p_rad=      atan(ts_Br/ts_Bp)
ts_p_deg=      180./pi*ts_p_rad
ts_etat_kmskpc= 1./3*ts_l_kpc*ts_v_kms
ts_Ralpha=      ts_alpk*ts_h/ts_etat
ts_Ralpha      [np.isnan(ts_Ralpha      )]=0.
ts_Ralphatot=   ts_alp*ts_h/ts_etat
ts_Ralphatot   [np.isnan(ts_Ralphatot   )]=0.
ts_ROmega=      ts_G*ts_h**2/ts_etat
ts_ROmega      [np.isnan(ts_ROmega      )]=0.
ts_D=           ts_Ralpha*ts_ROmega
ts_Dtot=        ts_Ralphatot*ts_ROmega
ts_Dc_anal=     -pi/32*(pi**2 +ts_Uz*ts_h/ts_etat)**2
ts_Dc_anal     [np.isnan(ts_Dc_anal     )]=0.
ts_DoverDc_anal=ts_D/ts_Dc_anal
ts_DoverDc_anal[np.isnan(ts_DoverDc_anal)]=0.
ts_DoverDtot   =ts_D/ts_Dtot
ts_DoverDtot   [np.isnan(ts_DoverDtot   )]=0.

#range_Gamma= np.zeros(nruns)
ts_Gamma=       np.zeros([nruns,ngmax,rows     ])
Gamma_mean=     np.zeros([nruns,ngmax          ])
ts_Gammatau=    np.zeros([nruns,ngmax,rows,cols])
ts_Gamma_Gyr=   np.zeros([nruns,ngmax,rows     ])
Gamma_mean_Gyr= np.zeros([nruns,ngmax          ])

dt=ts_t[0,0,1]-ts_t[0,0,0] #assumes equally spaced grid, same for all runs
for i in range(nruns):
  for j in range(int(ngals[i])):
    range_Gamma=[np.where((ts_t[i,j,:] > t1[i,j]) & (ts_t[i,j,:] < t2[i,j]))]
    ts_Gamma[i,j,:]= np.gradient(np.log(ts_Btot[i,j,:]),dt)
    Gamma_mean[i,j]=np.mean(ts_Gamma[i,j,range_Gamma])
    ts_Gammatau   [i,j,:]=Gamma_mean[i,j]*ts_tau[i,j,:]
    ts_Gamma_Gyr  [i,j,:]=ts_Gamma  [i,j,:]/t0_Gyr[i,j]
    Gamma_mean_Gyr[i,j  ]=Gamma_mean[i,j  ]/t0_Gyr[i,j]

#1) Bphi, Br vs r
plt.clf()
plt.xlim(0,r_disk_kpc[0,0])
plt.ylim(np.min(np.min([np.min(ts_Bp_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]),np.min(ts_Br_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])]) \
         -0.05*(np.max([np.max(ts_Bp_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]),np.max(ts_Br_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])]))), \
     np.max(1.05*(np.max([np.max(ts_Bp_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]),np.max(ts_Br_mkG[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])]))))
for i in range(nruns):
  for j in range(int(ngals[i])):
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_Br_mkG   [i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'r-' ,linewidth=i+1)
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_Bp_mkG   [i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'b--',linewidth=i+1)
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_Bzmod_mkG[i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'g-.',linewidth=i+1)
    plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')

plt.xlabel(r'$r\, [\mathrm{kpc}]$')
plt.ylabel(r'$B_\phi$, $B_r$, $B_z\, [\mathrm{mkG}]$')
plt.savefig('r_B_compare_multiple.eps', format='eps')

#2) alp_m, alp_k vs r
plt.clf()
plt.xlim(0,r_disk_kpc[0,0])
plt.ylim(np.min(np.min([np.min(ts_alpk_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]),np.min(ts_alpm_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])])\
    -0.05*(np.max([np.max(ts_alpk_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]),np.max(ts_alpm_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])]))), \
    np.max(1.05* (np.max([np.max(ts_alpk_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]]), \
    np.max(ts_alpm_kms[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])]))))
for i in range(nruns):
  for j in range(int(ngals[i])):
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_alpm_kms[i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'r-' ,linewidth=i+1)
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_alpk_kms[i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'b--',linewidth=i+1)
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],ts_alp_kms [i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],'g-.',linewidth=i+1)

plt.xlabel('$r$ [kpc]')
plt.ylabel(r'$\alpha_{\mathrm{m}}$, $\alpha_{\mathrm{k}}$, $\alpha\,[\mathrm{km\,s^{-1}}]$')
plt.savefig('alp_compare_multiple.eps', format='eps')

#3) -p vs r
plt.clf()
plt.xlim(0,r_disk_kpc[0,0])
plt.ylim([np.min(np.min(abs(ts_p_deg[:,:,nt,nxghost[0,0]+1:nx[0,0]-1-nxghost[0,0]])) \
                 -0.03*np.max(abs(ts_p_deg[:,:,nt,nxghost[0,0]+1:nx[0,0]-1-nxghost[0,0]]))), \
          np.max(1.03*np.max(abs(ts_p_deg[:,:,nt,nxghost[0,0]+1:nx[0,0]-1-nxghost[0,0]])))])
for i in range(nruns):
  for j in range(int(ngals[i])):
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-1-nxghost[i,j]],-ts_p_deg   [i,j,nt,nxghost[i,j]:nx[i,j]-1-nxghost[i,j]],'r-' ,linewidth=i+1)
    plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-1-nxghost[i,j]],-ts_pkin_deg[i,j,nt,nxghost[i,j]:nx[i,j]-1-nxghost[i,j]],'b--',linewidth=i+1)
                                                                                                        
plt.xlabel(r'$r\,[\mathrm{kpc}]$')
plt.ylabel(r'$-p$, $-p_{\mathrm{k}}\,(\mathrm{theor})\,[\mathrm{degrees}]$')
plt.savefig('p_compare_multiple.eps', format='eps')

#4) B vs t
plt.clf()

plt.subplot(211)
plt.xlim(0,10)
plt.ylim(np.min(1.1*np.min(np.log10(ts_Btot_mkG[:,:,0:n1[0,0]+1]))),np.max(1.1*np.max(np.log10(ts_Btot_mkG[:,:,0:n1[0,0]+1]))))
for i in range(nruns):
  for j in range(int(ngals[i])):
    plt.plot(ts_t_Gyr[i,j,0:n1[i,j]+1],np.log10(ts_Btot_mkG [i,j,0:n1[i,j]+1]),'r-',linewidth=i+1)
plt.plot([t1_Gyr,t1_Gyr],[-1.e5,1.e5],'k:')
plt.plot([t2_Gyr,t2_Gyr],[-1.e5,1.e5],'k:')

plt.subplot(212)
plt.xlim(0,10)
upprl=0.
for i in range(nruns):
  for j in range(int(ngals[i])):
    temp=1.1*np.max(ts_Gamma_Gyr[i,j,np.where((ts_t[i,j,:] > t1[i,j]) & (ts_t[i,j,:] < t2[i,j]))])
    if temp>upprl:
      upprl=temp
plt.ylim(0,upprl)
for i in range(nruns):
  for j in range(int(ngals[i])):
    plt.plot(ts_t_Gyr[i,j,0:n1[i,j]+1],ts_Gamma_Gyr[i,j,0:n1[i,j]+1],'r-',linewidth=i+1)

plt.ylabel(r'$(8\pi E_{\mathrm{tot}})^{1/2}\,[\mathrm{mkG}]$')
plt.plot([t1_Gyr,t1_Gyr],[-1.e5,1.e5],'k:')
plt.plot([t2_Gyr,t2_Gyr],[-1.e5,1.e5],'k:')
plt.xlabel(r'$t\,[\mathrm{Gyr}]$')
plt.ylabel(r'$\Gamma\,[\mathrm{Gyr}^{-1}]$')
plt.savefig('B_t_compare_multiple.eps', format='eps')

#5) profiles vs r
plt.clf()
fig=plt.figure(figsize=(20,10))
plt.subplots_adjust(left=0.04,right=0.97,bottom=0.08,top=0.93,wspace=0.3,hspace=0.5)

def prfplt(numr,numc,pos,arr,tit): #plots a profile of 'arr' vs r_kpc at position 'pos' and with title 'tit'
  for i in range(nruns):
    plt.subplot(numr,numc,pos)
    for j in range(int(ngals[i])):
      plt.plot(r_kpc[i,j,nxghost[i,j]:nx[i,j]-nxghost[i,j]],arr[i,j,nt,nxghost[i,j]:nx[i,j]-nxghost[i,j]],linewidth=i+1)
  plt.xlim(0,r_disk_kpc[0,0])
  plt.ylim(ymin=0)
  plt.ylim(ymax=np.max(1.5*np.max(arr[:,:,nt,nxghost[0,0]:nx[0,0]-nxghost[0,0]])))
  plt.xlabel(r'$r\,[\mathrm{kpc}]$')
  plt.title(tit)
  
numr=3 #number of rows in plot prfplt
numc=6 #number of cols in plot prfplt

prfplt(numr,numc, 1,ts_h_kpc       ,r'$h            \,[\mathrm{kpc}]$'                           ) #h(r)        disk scale height in kpc
prfplt(numr,numc, 2,ts_l_kpc       ,r'$l            \,[\mathrm{kpc}]$'                           ) #l(r)        turbulent outer scale in kpc
prfplt(numr,numc, 3,ts_v_kms       ,r'$v            \,[\mathrm{km\,s^{-1}}]$'                    ) #v(r)        rms turbulent velocity in km s^{-1}
prfplt(numr,numc, 4,ts_Uphi_kms    ,r'$U_{\phi}     \,[\mathrm{km\,s^{-1}}]$'                    ) #U_phi(r)    mean azimuthal velocity in km s^{-1}
prfplt(numr,numc, 5,ts_Uz_kms      ,r'$U_z          \,[\mathrm{km\,s^{-1}}]$'                    ) #U_z(r)      mean vertical velocity in km s^{-1}
prfplt(numr,numc, 6,ts_tau_Myr     ,r'$\tau=l/v     \,[\mathrm{Myr}]$'                           ) #tau(r)      turb correl time l/v in Myr^{-1}
prfplt(numr,numc, 7,ts_H           ,r'$H\equiv h/l$'                                             ) #H(r)        Dim/less disk scale height H=h/l
prfplt(numr,numc, 8,ts_Phi         ,r'$\Phi\equiv\Omega\tau$'                                    ) #Phi(r)      Inv Rossby=0.5*Coriolis, Phi=Om*tau
prfplt(numr,numc, 9,ts_V           ,r'$V\equiv U_z/v$'                                           ) #V(r)        Dim/less vertical mean veloc V=U_z/v
prfplt(numr,numc,10,ts_n_cm3       ,r'$n            \,[\mathrm{cm^{-3}}]$'                       ) #n(r)        gas number density in cm^{-3}
prfplt(numr,numc,11,ts_Beq_mkG     ,r'$B_\mathrm{eq}  \,[\mathrm{m\mu G}]$'                      ) #B_eq(r)     equipart mag fld strength in microG
prfplt(numr,numc,12,-ts_G_Gyr      ,r'$-S\equiv r\partial\Omega/\partial r\,[\mathrm{Gyr^{-1}}]$') #G(r)        radial shear r*d\Om/dr in Gyr^{-1}
prfplt(numr,numc,13,ts_etat_cm2s   ,r'$\eta_\mathrm{t}\,[\mathrm{10^{26}cm^2s^{-1}}]$'           ) #etat(r)     turb diffusivity in cm^2 s^{-1}
prfplt(numr,numc,14,-ts_D          ,r'$-D$'                                                      ) #D(r)        Dynamo num D=alpk*G*h^3/etat^2
prfplt(numr,numc,15,-ts_Dc_anal    ,r'$-D_\mathrm{c}}\mathrm{(analytical)}$'                     ) #D_c(r)      Analytical estimate of crit Dyn num
prfplt(numr,numc,15,-ts_Dtot       ,r'$-D_\mathrm{c}\mathrm{(analytical)},\,D_{\mathrm{tot}}$'   ) #D_tot(r)    Effective total dyn number
prfplt(numr,numc,16,ts_DoverDtot   ,r'$D/D_\mathrm{c}(D/D_{\mathrm{tot}}$'                       ) #D(r)/D_tot(r) 
prfplt(numr,numc,16,ts_DoverDc_anal,r'$D/D_\mathrm{c}(\mathrm{analytical}),\,D/D_{\mathrm{tot}}$') #D(r)/D_c(r) 
prfplt(numr,numc,17,ts_Ralpha      ,r'$R_{\alpha}$'                                              ) #R_alpha     Dim/less alp_k R_alp=alpk*h/etat
prfplt(numr,numc,17,ts_Ralphatot   ,r'$R_{\alpha},\,R_{\alpha,\mathrm{tot}}$'                    ) #R_alpha,tot Dim/less alpha R_alp,tot=alp*h/etat
prfplt(numr,numc,18,-ts_ROmega     ,r'$-R_{\Omega}$'                                             ) #R_Omega     Dim/less shear R_Om=S*h^2/etat

plt.savefig('profiles_r_compare_multiple.eps', format='eps')

os.chdir(programdirec)
print 'current working directory=',os.getcwd()
