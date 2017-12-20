!I.A. Grant Wilson
!05/Jan/2017
!DE_full_grid_normal.f90

!MIT License

!Copyright (c) [2017] [Grant Wilson - ]

!Permission is hereby granted, free of charge, to any person obtaining a copy
!of this software and associated documentation files (the "Software"), to deal
!in the Software without restriction, including without limitation the rights
!to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!copies of the Software, and to permit persons to whom the Software is
!furnished to do so, subject to the following conditions:

!The above copyright notice and this permission notice shall be included in all
!copies or substantial portions of the Software.

!THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
!IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
!FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
!AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
!LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
!OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
!SOFTWARE.

!Maximising the revenue available to a storage system that can buy from the spot price market.
!Code adapted from paper: DOI 10.1039/C2EE02419E Edward Barbour, Grant Wilson et al.
!make sure power is in MW and prices are per MWh! This is different from all previous versions
!which were in kW and kWh

!! e_to_store is within the store****!!!
!! output_to_grid is outwith the store (requires efficiencies)***!!
!! grid_constraint has efficiencies - so must be outwith the store
!! ELI and ELO do not have efficiencies - so must be within the store


PROGRAM  Full_grid_1

IMPLICIT NONE

INTEGER, PARAMETER :: length = 61088
INTEGER, PARAMETER :: num_rev_output = 10000

CHARACTER::market*2
INTEGER :: i, n, numtrials, step_size, period_1, period_2
INTEGER :: temp, accept, ifill, window, year, scenario, Cap_max, num_rev_count
INTEGER :: R, i_temp, i_temp2, accept_count, revenue_iterations(num_rev_output)
INTEGER :: flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8, flags_total
INTEGER :: flag1_tot, flag2_tot, flag3_tot, flag4_tot, flag5_tot, flag6_tot, flag7_tot, flag8_tot

REAL (kind = 8) :: prices(length)
REAL (kind = 8) :: time(length)
REAL (kind = 8) :: t(length)
REAL (kind = 8) :: timecol(length)
REAL (kind = 8) :: Old_OTG(length), Old_ETS(length), Output_to_grid(length), grid_constraint(length), Cumsum_Rev(length)
REAL (kind = 8) :: old_test_energy(length), test_energy(length), new_test_energy(length)
REAL (kind = 8) :: MWh_price(length), Revenue(length), revenue_iterations_value(num_rev_output)
REAL (kind = 8) :: E_stored(length), E_to_store(length)
!REAL (kind = 8) :: Cap_max
REAL (kind = 8) :: HARVEST
REAL (kind = 8) :: prob_accept
REAL (kind = 8) :: dx, dx1, dx2, dx3, dx4, dx5
REAL (kind = 8) :: eta, eta_in, eta_out
REAL (kind = 8) :: tau
REAL (kind = 8) :: PLI, PLO, ELI, ELO
REAL (kind = 8) :: nn_bias
REAL (kind = 8) :: limit
REAL (kind = 8) :: time_loss, time_loss_2, time_loss_3, time_loss_4
REAL (kind = 8) :: old_Rev, new_Rev, New_quantity, Original_quantity
REAL (kind = 8) :: Rev1, Rev2, new_OtG_1, new_OtG_2
REAL (kind = 8) :: m, m2, time1, time2, total_cpu_time

call cpu_time ( time1 )
READ (*,*) eta_in
READ (*,*) eta_out
READ (*,*) Cap_max
READ (*,*) PLI
READ (*,*) PLO
READ (*,*) tau
READ (*,*) numtrials
READ (*,*) step_size
READ (*,*) nn_bias
READ (*,*) market
READ (*,*) year


ELI = PLI
ELO = -PLO
!ELI = ELI*eta_in
!ELO = ELO*eta_out
IF (step_size>ELI) THEN
step_size = ELI
END IF

!****time file has to be in the correct place on the server***
OPEN(UNIT=22, FILE = '../../time_input_files/2050_time_period_DE_GB_hour.txt')
OPEN(UNIT=25, FILE = 'price.txt')
OPEN(UNIT=77, FILE = '../../time_input_files/2050_time_DE_GB_hour.txt')

!RECL increases width of output file before Linefeed
!OPEN(UNIT=73, FILE = 'scenario_moves.txt', RECL = 550)
!OPEN(UNIT=74, FILE = 'all_dx.txt', RECL = 750) 
OPEN(UNIT=76, FILE = 'revenue_iterations.txt', RECL = 550) 

READ(22,*) time
READ(25,*) prices
READ(77,*) timecol

CLOSE(22)
CLOSE(25)
CLOSE(77)


old_Rev = SUM(Revenue)
!change window to less than a year (17520) e.g. first month only is 48 * 31 = 1488
!first week is 336 and first 4 days are 192
window=61088 
accept_count=1
num_rev_count=1
flag1_tot=0; flag2_tot=0; flag3_tot=0; flag4_tot=0; flag5_tot=0; flag6_tot=0; flag7_tot=0; flag8_tot=0
flags_total=0

!Initialise arrays
DO i = 1,length
t(i) = time(i)
Cumsum_Rev(i)=0
old_test_energy(i) = 0
MWh_price(i) = prices(i) !Prices are in Eur/MWh or GBP/MWh depending on market

grid_constraint(i) = ELI   !set as same limit as energy limit in
!Initialise other arrays
E_stored(i) = 0
E_to_store(i) = 0
Revenue(i) = 0
Output_to_grid(i) = 0
Old_OTG(i)=0   !for debugging purposes
Old_ETS(i)=0   !for debugging purposes
END DO

DO i = 1, num_rev_output
revenue_iterations(i)=0
revenue_iterations_value(i)=0
END DO

!Start the Monte carlo trials...

CALL RANDOM_SEED  !this sets the pseudo-random number seed to a different starting point
DO n = 1,numtrials

!Set up the variable that accepts or declines proposed changes
accept = 0 !Default is to decline

DO
!Choose a random amount to shift, either +ve or -ve.
CALL RANDOM_NUMBER(HARVEST)
dx = (step_size - ((HARVEST) * step_size * 2))
!Select a random period
CALL RANDOM_NUMBER(HARVEST)
period_1 = NINT(1 + ((HARVEST) * (window-1)))
!Biasing towards nearest neighbours and making sure period_1 /= period_2

CALL RANDOM_NUMBER(HARVEST)
!period_2 = NINT(1 + ((HARVEST) * (window-1)))  !previous probability method - takes a lot of cpu time
!period_2 = period_1 - (nn_bias - ((HARVEST) * nn_bias * 2)) ! test check - this did not work with big nn_bias
period_2 = NINT(period_1 + ((HARVEST) * nn_bias))
accept = 1

!make sure period_1 is smaller
!IF (period_2<period_1) THEN
!temp = period_1
!period_1 = period_2
!period_2 = temp
!END IF

!IF (period_2 < (period_1 + nn_bias)) THEN; accept = 1; END IF
IF (period_1 == period_2) THEN; accept=0; ENDIF !make sure periods are not equal
!IF (period_2 < 1) THEN; accept=0; ENDIF
IF (period_2 > (window)) THEN; accept=0; ENDIF 
!Check for the exit command
IF (accept==1) EXIT 

END DO

!set the default back to decline
accept = 0
flag1=0; flag2=0; flag3=0; flag4=0; flag5=0; flag6=0; flag7=0; flag8=0




time_loss = EXP((t(period_1)-t(period_2))/tau)


IF (dx>0) THEN


!scenario 1
            IF(E_to_store(period_1)>=0 .AND. E_to_store(period_2)>0) THEN
                scenario=1
                dx1 = (ELI-E_to_store(period_1))                                               !storage limit in
                dx2 = E_to_store(period_2)/time_loss                                           !storage limit out
                dx3 = (grid_constraint(period_1)+Output_to_grid(period_1))*eta_in              !grid limit in
                dx4 = -Output_to_grid(period_2)/(eta_out*time_loss)                            !grid limit out
                dx = MIN (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 2
            IF(E_to_store(period_1)>=0 .AND. E_to_store(period_2)<=0) THEN
                scenario=2
                dx1 = (ELI-E_to_store(period_1))                                               !storage limit in
                dx2 = (E_to_store(period_2)-ELO)/time_loss                                     !storage limit out
                dx3 = (grid_constraint(period_1)+Output_to_grid(period_1))*eta_in              !grid limit in
                dx4 = (grid_constraint(period_2)-Output_to_grid(period_2))/(eta_out*time_loss) !grid limit out
                dx = MIN (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 3
            IF(E_to_store(period_1)<0 .AND. E_to_store(period_2)>0) THEN
                scenario=3
                dx1 = -E_to_store(period_1)                                                    !storage limit in
                dx2 = E_to_store(period_2)/time_loss                                           !storage limit out
                dx3 = Output_to_grid(period_1)*eta_in                                          !grid limit in
                dx4 = -Output_to_grid(period_2)/(eta_out*time_loss)                            !grid limit out
                dx = MIN (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 4
            IF(E_to_store(period_1)<0 .AND. E_to_store(period_2)<=0) THEN
                scenario=4
                dx1 = -E_to_store(period_1)                                                     !storage limit in
                dx2 = (E_to_store(period_2)-ELO)/time_loss                                      !storage limit out
                dx3 = Output_to_grid(period_1)*eta_in                                           !grid limit in
                dx4 = (grid_constraint(period_2)-Output_to_grid(period_2))/(eta_out*time_loss)  !grid limit out
                dx = MIN (dx, dx1, dx2, dx3, dx4)
            END IF

!This section checks the available space for moving energy between periods 1 and 2
            i_temp = MAXLOC(E_stored(period_1:(period_2-1)), DIM=1)
            i_temp2 = period_1+i_temp-1
                    
            time_loss_3 = EXP((time(period_1)-time(i_temp2))/tau)
            m2 = Cap_max - MAXVAL(E_stored(period_1:(period_2 - 1)))
            dx5 = (m2/time_loss_3)
                                
            IF(dx5<0) THEN
            dx5 = 0
            END IF
            
            dx = MIN(dx, dx5)
            IF (dx < 0) THEN; dx=0; ENDIF


!This section checks the revenue of the potential move - accept if revenue increases
!Works out the new values Output_to_grid would be if the move is accepted

accept = 0        
            IF(scenario==1) THEN
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in
            END IF

            IF(scenario==2) THEN
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out)
            END IF

            IF(scenario==3) THEN
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in
            END IF

            IF(scenario==4) THEN
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out)
            END IF

Rev1=Output_to_grid(period_1)*MWh_price(period_1) + Output_to_grid(period_2)*MWh_price(period_2)
Rev2=new_OtG_1*MWh_price(period_1) + new_OtG_2*MWh_price(period_2)

IF(new_OtG_1<=grid_constraint(period_1) .AND. new_OtG_2>=-grid_constraint(period_2)) THEN
IF(Rev2 > Rev1) THEN; accept=1; END IF
END IF

ELSE 

!dx<0



!scenario 5
            IF(E_to_store(period_1)>0 .AND. E_to_store(period_2)>=0) THEN
                scenario=5
                dx1 = -E_to_store(period_1)                                                      !storage limit in
                dx2 = -(ELI-E_to_store(period_2))/time_loss                                      !storage limit out
                dx3 = Output_to_grid(period_1)*eta_in                                            !grid limit in
                dx4 = -((grid_constraint(period_2)+Output_to_grid(period_2))/(eta_out*time_loss))  !grid limit out
                dx = MAX (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 6
            IF(E_to_store(period_1)>0 .AND. E_to_store(period_2)<0) THEN
                scenario=6
                dx1 = -E_to_store(period_1)                                                     !storage limit in
                dx2 = E_to_store(period_2)/time_loss                                      !storage limit out
                dx3 = Output_to_grid(period_1)*eta_in                                          !grid limit in
                dx4 = -(Output_to_grid(period_2))/(eta_out*time_loss)  !grid limit out
                dx = MAX (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 7
            IF(E_to_store(period_1)<=0 .AND. E_to_store(period_2)>=0) THEN
                scenario=7
                dx1 = (ELO-E_to_store(period_1))                                                     !storage limit in
                dx2 = -(ELI-E_to_store(period_2))/time_loss                                     !storage limit out
                dx3 = -(grid_constraint(period_1)-Output_to_grid(period_1))*eta_in                                          !grid limit in
                dx4 = -((grid_constraint(period_2)+Output_to_grid(period_2))/(eta_out*time_loss))  !grid limit out
                dx = MAX (dx, dx1, dx2, dx3, dx4)
            END IF
            
!scenario 8
            IF(E_to_store(period_1)<=0 .AND. E_to_store(period_2)<0) THEN
                scenario=8
                dx1 = (ELO-E_to_store(period_1))                                                     !storage limit in
                dx2 =  E_to_store(period_2)/time_loss                                            !storage limit out
                dx3 = -(grid_constraint(period_1)-Output_to_grid(period_1))*eta_in                                          !grid limit in
                dx4 = -(Output_to_grid(period_2))/(eta_out*time_loss)  !grid limit out  !grid limit out
                dx = MAX (dx, dx1, dx2, dx3, dx4)
            END IF

!This section checks the available space for moving energy between periods 1 and 2
            i_temp = MINLOC(E_stored(period_1:(period_2-1)), DIM=1)
            i_temp2 = period_1+i_temp-1
                    
            time_loss_4 = EXP((time(period_1)-time(i_temp2))/tau)
            m = MINVAL(E_stored(period_1:(period_2 - 1)))
            dx5 = (-m/time_loss_4)!+(1e-8)
                                
            IF(dx5>0) THEN
            dx5 = 0
            END IF
            
            dx = MAX(dx, dx5)

            IF (dx > 0) THEN; dx=0; ENDIF


!This section checks the revenue of the potential move - accept if revenue increases
!Works out the new values Output_to_grid would be if the move is accepted

accept = 0
            IF(scenario==5) THEN 
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in
            END IF
           
            IF(scenario==6) THEN 
                new_OtG_1=Output_to_grid(period_1)-dx/eta_in
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out)
            END IF
 
            IF(scenario==7) THEN 
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss)/eta_in
            END IF

            IF(scenario==8) THEN 
                new_OtG_1=Output_to_grid(period_1)-dx*eta_out
                new_OtG_2=Output_to_grid(period_2)+(dx*time_loss*eta_out)
            END IF

Rev1=Output_to_grid(period_1)*MWh_price(period_1) + Output_to_grid(period_2)*MWh_price(period_2)
Rev2=new_OtG_1*MWh_price(period_1) + new_OtG_2*MWh_price(period_2)
IF(new_OtG_1<=grid_constraint(period_1) .AND. new_OtG_2>=-grid_constraint(period_2)) THEN 
IF(Rev2 > Rev1) THEN; accept=1; END IF 
END IF 


END IF  

IF(dx==0) THEN; accept = 0; END IF 


!now if the move is accepted then all the variables can be updated..
IF(accept == 1) THEN 
        
        time_loss_2 = EXP((time(period_1)-time(period_2))/tau)

!update energy in store at periods 1 and 2       
        E_to_store(period_1) = E_to_store(period_1) + dx
        E_to_store(period_2) = E_to_store(period_2) - dx*time_loss_2
!update output to grid at periods 1 and 2
        Output_to_grid(period_1)=new_OtG_1
        Output_to_grid(period_2)=new_OtG_2

!update state of charge between periods 1 and 2
DO ifill = period_1,(period_2 - 1) 
E_stored(ifill)=E_stored(ifill)+(dx*EXP((time(period_1)-time(ifill))/tau)) 
END DO 
 

!Check for flags in debugging
IF(scenario==1) THEN; IF (dx < 0) THEN; flag1=1; flag1_tot=flag1_tot+1; ENDIF; END IF
IF(scenario==2) THEN; IF (dx < 0) THEN; flag2=1; flag2_tot=flag2_tot+1; ENDIF; END IF 
IF(scenario==3) THEN; IF (dx < 0) THEN; flag3=1; flag3_tot=flag3_tot+1; ENDIF; END IF
IF(scenario==4) THEN; IF (dx < 0) THEN; flag4=1; flag4_tot=flag4_tot+1; ENDIF; END IF
IF(scenario==5) THEN; IF (dx > 0) THEN; flag5=1; flag5_tot=flag5_tot+1; ENDIF; END IF
IF(scenario==6) THEN; IF (dx > 0) THEN; flag6=1; flag6_tot=flag6_tot+1; ENDIF; END IF
IF(scenario==7) THEN; IF (dx > 0) THEN; flag7=1; flag7_tot=flag7_tot+1; ENDIF; END IF
IF(scenario==8) THEN; IF (dx > 0) THEN; flag8=1; flag8_tot=flag8_tot+1; ENDIF; END IF


END IF

flags_total=(flag1_tot+flag2_tot+flag3_tot+flag4_tot+flag5_tot+flag6_tot+flag7_tot+flag8_tot)


!The following sections allow output of data for debugging etc.
IF (accept==1) THEN
!write to files the first accepted 2000 values for total revenue, then every 100th - used for figures to show trending to final value
IF (accept_count<2000) THEN; revenue_iterations(num_rev_count) =  n;  revenue_iterations_value(num_rev_count) = SUM(Output_to_grid*MWh_price); num_rev_count=num_rev_count+1; ENDIF
IF ((MOD((accept_count-2000),100)==1) .AND. (num_rev_count<num_rev_output)) THEN; revenue_iterations(num_rev_count) =  n;  revenue_iterations_value(num_rev_count) = SUM(Output_to_grid*MWh_price); num_rev_count=num_rev_count+1; ENDIF

!next section outputs to all_dx.txt file for debugging
!IF (accept_count==1) THEN;
!WRITE(74,*) ' n ', 'accept', ' period_1 ', ' period_2 ', ' scenario ', ' dx ', ' dx1 ', ' dx2 ', ' dx3 ', ' dx4 ', ' dx5 ', ' TL1 ', ' TL2 ', ' TL4 ', ' Old_OTG1 ', ' Old_OTG2 ', ' OTGP1 ', ' OTGP2 ', ' new_OtG_1 ', ' new_OtG_2 ', ' Old_ETS1 ', ' Old_ETS2 ', ' ETS1 ',  ' ETS2 ', ' Rev2-Rev1 ', ' gc1 ', ' gc2 ',  ' flag1 ', ' flag2 ', ' flag3 ', ' flag4 ', ' flag5 ', ' flag6 ', ' flag7 ', ' flag8 ', ' flags_total ' !all_dx.txt file see line ~68
!ENDIF
!WRITE(74,*),  n, accept, period_1, period_2, scenario, dx, dx1, dx2, dx3, dx4, dx5, time_loss, time_loss_2, time_loss_4, Old_OTG(period_1), Old_OTG(period_2), Output_to_grid(period_1), Output_to_grid(period_2), new_OtG_1, new_OtG_2, Old_ETS(period_1), Old_ETS(period_2),  E_to_store(period_1),  E_to_store(period_2), Rev2-Rev1, grid_constraint(period_1), grid_constraint(period_2), flag1, flag2, flag3, flag4, flag5, flag6, flag7, flag8, flags_total  !all_dx.txt file see line ~68
!140 format (5a10, 9a5, 22a11)
!150 format (5i15,22f15.3,9i6)
!DO i = 1,window
!use line below to make output data used to prepare video
!WRITE (73,*), n, accept_count, scenario, E_stored(i), Rev2-Rev1
!END DO


accept_count = accept_count+1
ENDIF

Old_OTG(period_1)=Output_to_grid(period_1)
Old_OTG(period_2)=Output_to_grid(period_2)
Old_ETS(period_1)=E_to_store(period_1)
Old_ETS(period_2)=E_to_store(period_2)

 
END DO ! end of main program loop





!The following sections output data for graphing and debugging

OPEN(UNIT=60, FILE = 'Run_deets.txt')
OPEN(UNIT=72, FILE = 'output_for_figures.csv', RECL = 550)
OPEN(UNIT=78, FILE = 'Revenue.csv')

DO i = 1,window
Revenue(i) = Output_to_grid(i)*MWh_price(i)
IF (i==1) THEN; Cumsum_Rev(i)=Revenue(i); ENDIF
IF (i>1) THEN; Cumsum_Rev(i)=Cumsum_Rev(i-1)+Revenue(i); ENDIF
IF (i==1) THEN; WRITE(72,FMT='(a2,a11,a2,a8,a2,a11,a2,a17,a2,a15,a2,a10,a2,a5,a2,a8)'),market,'_Date_Time,',market,'_period,',market,'_MWh_price,',market,'_E_to_store(MWh),',market,'_E_stored(MWh),',market,'_OTG(MWh),',market,'_Rev,',market,'_Cum_Rev'; ENDIF
WRITE(72, FMT='(f12.6,a,i,a,f10.2,a,f14.5,a,f14.5,a,f14.5,a,f14.2,a,f14.2)'),timecol(i),',',i,',',MWh_price(i),',',E_to_store(i),',',E_stored(i),',',Output_to_grid(i),',',Revenue(i),',',Cumsum_Rev(i)
IF (i==1) THEN; WRITE(78,FMT='(a2,a1,f4.2,a1,f5.0,a1,f5.0,a1,I6.6)'),market,'_',eta_in*eta_out,'_',PLI,'_',PLO,'_',Cap_max; ENDIF !Cap_max changed from floating point to Integer to allow leading zeros
WRITE(78, FMT='(f14.2)'), Revenue(i)
END DO
CLOSE(72);CLOSE(78);!CLOSE(73)


new_Rev = SUM(Revenue)
Original_quantity = SUM(new_test_energy)
New_quantity = SUM(Output_to_grid)


OPEN(UNIT=61, FILE = 'upper_boundary.txt')
WRITE(61, *), new_Rev
CLOSE(61)



OPEN(UNIT=81, FILE = 'flag_scenarios.txt', RECL = 550)
WRITE(81, FMT='(9a12)'), 'flag1', 'flag2', 'flag3', 'flag4', 'flag5', 'flag6', 'flag7', 'flag8', 'flags_total'
WRITE(81, FMT='(9i12)'), flag1_tot, flag2_tot, flag3_tot, flag4_tot, flag5_tot, flag6_tot, flag7_tot, flag8_tot, flags_total
CLOSE(81)

!OPEN(UNIT=82, FILE = 'flags_total.txt')
!WRITE(82, *), flags_total
!CLOSE(82)

DO i = 1,num_rev_output
IF (revenue_iterations(i) > 0) THEN
WRITE (76, FMT='(i,f12.2)') revenue_iterations(i),  revenue_iterations_value(i)
END IF
END DO
CLOSE(76)

call cpu_time ( time2 )

110 format (a45,f15.4)
120 format (a45,i15)
130 format (a45,a15)
WRITE (60,110) '(time independent) round trip efficiency = ', eta_in*eta_out
WRITE (60,110) '(time independent) efficiency in (1=100%) = ', eta_in
WRITE (60,110) '(time independent) efficiency out (1=100%) = ', eta_out
WRITE (60,120) 'max capacity (MWh) = ', Cap_max
WRITE (60,110) 'Power Limit In (MW) = ', PLI
WRITE (60,110) 'Power Limit Out (MW) = ', PLO
IF (tau < 10000) THEN; WRITE (60,110) 'storage time constant (hours)  = ', tau; END IF
IF (tau >= 10000) THEN; WRITE (60,FMT='(a45,es15.3e2)') 'storage time constant (hours)  = ', tau; END IF
WRITE (60,120) 'number of trials = ', numtrials
WRITE (60,120) 'step_size = ', step_size
WRITE (60,110) 'nearest neighbour Bias = ', nn_bias
WRITE (60,130) 'market = ', market
WRITE (60,120) 'year = ', year
!140 and 150 used previously
WRITE (60,110) 'Old Revenue = ', old_Rev
WRITE (60,110) 'new Revenue = ', new_Rev
WRITE (60,110) 'Original energy output = ', Original_quantity
WRITE (60,110) 'New energy output = ', New_quantity
WRITE (60,120) 'Flags total = ', flags_total
WRITE (60,110) 'total cpu time = ', time2-time1

CLOSE (60)

!STOP

END PROGRAM Full_grid_1

