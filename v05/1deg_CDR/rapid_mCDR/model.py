import numpy as np
import xarray as xr

def TDMASolve(a, b, c, d):
    """
    solver for a tridiagonal matrix.
    [a b c] x = d
    a, b, c are the diagonals
    d is the right hand side
    """


    n = len(a)
    ac, bc, cc, dc = map(np.array, (a, b, c, d))
    xc = []
    for j in range(1, n):
        if(bc[j - 1] == 0):
            ier = 1
            return
        ac[j] = ac[j]/bc[j-1]
        bc[j] = bc[j] - ac[j]*cc[j-1]
    if(b[n-1] == 0):
        ier = 1
        return
    for j in range(1, n):
        dc[j] = dc[j] - ac[j]*dc[j-1]
    dc[n-1] = dc[n-1]/bc[n-1]
    for j in range(n-2, -1, -1):
        dc[j] = (dc[j] - cc[j]*dc[j+1])/bc[j]
    return dc



def solve_prog_eqn(var,w,k,z,zl,dt,surf_forc,central_diff=False):
    """
    implicit solution of DIC/ALK equation
    """
    a=np.zeros(var.shape)
    b=np.zeros(var.shape)
    c=np.zeros(var.shape)
    
    d=var
    d[0]=d[0]+surf_forc*dt/(zl[1]-zl[0])


    A=dt/(zl[2:-1]-zl[1:-2])
    f1=A*w[1:-2]
    f2=A*w[2:-1]

    # factors for upwind approximation (alpha_k,alpha_kp1=0.5 --> central differencing)    
    alpha_k=np.heaviside(w[1:-2],0)
    alpha_kp1=np.heaviside(w[2:-1],0)

    f3=k[1:-2]*dt/(z[1:-1]-z[0:-2])/(zl[2:-1]-zl[1:-2])
    f4=k[2:-1]*dt/(z[2:]-z[1:-1])/(zl[2:-1]-zl[1:-2])

    a[1:-1]=-f1*alpha_k-f3
    b[1:-1]=1+f2*alpha_kp1-f1*(1-alpha_k)+f3+f4
    c[1:-1]=f2*(1-alpha_kp1)-f4

    # boundary conditions
    # zero tendency at the bottom
    b[-1]=1
     
    # zero flux at the surface

    alpha_0=np.heaviside(w[1],0)

    f1s=dt/(zl[1]-zl[0])*w[1]
    f2s=k[1]*dt/(z[1]-z[0])/(zl[1]-zl[0])
    b[0]=1+f1s*alpha_0+f2s
    c[0]=f1s*(1-alpha_0)-f2s

    var_new=TDMASolve(a, b, c, d)
    return var_new



def rapid_mcdr_model(dta):
    """
    integrate rapid-mcdr prognostic equations for dAlk and dDIC in time 
    (the time step is defined in the input data)

    input: dta .... forcing_data
    xarray.dataset with the following variables:
       time (time) ... the model is tested with a daily time-step 
       Zl (k_l) ... vertical coordinate of lower interface (m, negative values) 
       w (time,k_l) ... vertical velocity in m/s
       k_diff (time,k_l) ... diffusivity for DIC/ALK in m^2/s
       k_surf (time) ... piston velocity
       dpco2_ov_dalk (time) ... dpCO2/dALK|DIC in atm /(mol ALK/kg)
       dpco2_ov_ddic (time) .... dpCO2/dDIC|ALK in atm /(mol C/kg)
       alk_forcing (time) ... alkalinity forcing in meq s-1
       area (time) ... area of forcing in m^2
       siarea (time) .... ice covered area [0-1]

    output: raid_mcdr_results
    xarray, see metadata for info   

    """

    # 1s (useful to convert data to seconds)
    dts=np.timedelta64(1,'s')

    #if spacial_averaging_forcing=='pco2':
    #    dta=xr.open_dataset(path_pert+'/1d_forcing/model_forcing_data_daily_pco2.nc')
    #elif spacial_averaging_forcing=='default':
    #    dta=xr.open_dataset(path_pert+'/1d_forcing/model_forcing_data_daily.nc')
    #else:
    #    raise Exception('Unknown spacial_averaging_forcing')    
    
    time=dta.time
    n_time=len(dta.time.values)

    # z_l  depth of the mid layers - first value will be zero 
    zl=np.abs(dta.Zl.values)
    # z - depth of full layers
    z=0.5*(zl[1:]+zl[0:-1])
    dz=zl[1:]-zl[0:-1]
    n_z=len(z)
    # vertical velocity and diffusivity on time*z_l
    w=-dta.w.values
    k_diff=dta.k_diff.values

    # all units of d_dic and d_alk are in uM = meq m-3 
    # 1 uM = 10^-3 mol m^-3
    d_dic=np.zeros((n_time,n_z))*np.nan
    d_alk=np.zeros((n_time,n_z))*np.nan

    f_co2=np.zeros(n_time)*np.nan
    
    d_dic[0,:]=0
    d_alk[0,:]=0
    f_co2[0]=0

    for t in range(n_time-1):
    # integrate from time t --> t+1
        # dt in seconds
        dt=((time[t+1]-time[t])/dts).values

        forc_co2=-(dta['dpco2_ov_ddic'].values[t]*d_dic[t,0]+dta['dpco2_ov_dalk'].values[t]*d_alk[t,0])*dta['k_surf'].values[t]*(1-dta['siarea'][t])
        
        # alk_forcing is the total forcing - need to divide by area to get it in the per-area units
        forc_alk=dta['alk_forcing'].values[t]/dta['area'].values[t]

        d_dic[t+1,:]=solve_prog_eqn(d_dic[t,:],w[t,:],k_diff[t,:],z,zl,dt,forc_co2) #,central_diff=central_diff)
        d_alk[t+1,:]=solve_prog_eqn(d_alk[t,:],w[t,:],k_diff[t,:],z,zl,dt,forc_alk) #,central_diff=central_diff)

        f_co2[t]=forc_co2

    res = xr.Dataset(coords={'time':('time',dta.time.values),'z':('z',z)})
    res['f_co2']=('time',f_co2)
    res['f_co2']=res['f_co2'].assign_attrs({'long_name':'surface CO2 flux','units':'?'})
    
    res['dic']=(('time','z'),d_dic)
    res['dic']=res['dic'].assign_attrs({'long_name':'dDIC','units':'uM'})
   
    res['alk']=(('time','z'),d_alk)
    res['alk']=res['alk'].assign_attrs({'long_name':'dALK','units':'uM'})
    
    # spaceally-integrated delta DIC/ALK
    res['dDIC_profile']=res['dic']*dta['area'] 
    res['dALK_profile']=res['alk']*dta['area']


    cdr_potential=dta['dpco2_ov_dalk']/dta['dpco2_ov_ddic']
    res['cdr_potential']=cdr_potential

    # that does make sense only for continious experiments
    #res['X']=-f_co2/dta['alk_forcing']*dta['area']/cdr_potential
    #res['X']=res['X'].assign_attrs({'long_name':'Surface CO2 flux','units':'normalized'})
    
    res['dpco2_ov_ddic']=dta['dpco2_ov_ddic']
    res['dpco2_ov_dalk']=dta['dpco2_ov_dalk']
    res['k_surf']=dta['k_surf']
    res['area']=dta['area']
    res['alk_forcing']=dta['alk_forcing']
    
    res['dz']=('z',dz)
    res['dz']=res['dz'].assign_attrs({'long_name':'delta z','units':'m'})

    return res  