MODULE eos_newt_functions

implicit none

CONTAINS

    FUNCTION newt_eos(x)
        use nrtype
        use Driver_interface, ONLY: Driver_abortFlash
        use Eos_interface, ONLY: Eos
        use Eos_helmData, ONLY: nr_eos_mode, nr_eos_eint, nr_eos_pres, nr_eos_entr, nr_eos_deriv, nr_eos_k, &
            eos_smallt, eos_larget, nr_eos_scale, nr_eos_limited, nr_eos_lim_cnt
        use eos_vecData, ONLY: denRow, tempRow, etotRow, ptotRow, stotRow, dptRow, detRow, dpdRow, dstRow, dsdRow
        use Grid_data, ONLY: gr_smallrho
        implicit none
#include "constants.h"
#include "Eos.h"

        real, dimension(:), intent(in) :: x
        real, dimension(size(x)) :: newt_eos
        double precision :: tlim, dlim, sclim
        nr_eos_limited = .false.

        !print *, eos_smallt, eos_larget
        if (nr_eos_mode .eq. MODE_DENS_EI .or. nr_eos_mode .eq. MODE_DENS_PRES) then
            tlim = x(1)*nr_eos_scale(1)
            if (tlim .lt. eos_smallt) then
                tlim = eos_smallt
                nr_eos_limited = .true.
            endif
            if (tlim .gt. eos_larget) then
                tlim = eos_larget
                nr_eos_limited = .true.
            endif
            tempRow(nr_eos_k) = tlim
        else
            dlim = x(1)*nr_eos_scale(1)
            if (dlim .lt. gr_smallrho) then
                dlim = gr_smallrho
                nr_eos_limited = .true.
            endif 
            denRow(nr_eos_k) = dlim
        endif
        if (nr_eos_mode .eq. MODE_PRES_ENTR) then
            if (size(x) .ne. 2) then
                call Driver_abortFlash('x wrong size!')
            endif
            tlim = x(2)*nr_eos_scale(2)
            if (tlim .lt. eos_smallt) then
                tlim = eos_smallt
                nr_eos_limited = .true.
            endif
            if (tlim .gt. eos_larget) then
                tlim = eos_larget
                nr_eos_limited = .true.
            endif
            tempRow(nr_eos_k) = tlim
        endif
        call eos_helm(nr_eos_k, nr_eos_k)
        if (nr_eos_mode .eq. MODE_DENS_PRES) then
            newt_eos(1) = (ptotRow(nr_eos_k) - nr_eos_pres)/abs(nr_eos_pres)
            nr_eos_deriv(1,1) = dptRow(nr_eos_k)/abs(nr_eos_pres)*nr_eos_scale(1)
        elseif (nr_eos_mode .eq. MODE_DENS_EI) then
            newt_eos(1) = (etotRow(nr_eos_k) - nr_eos_eint)/abs(nr_eos_eint)
            nr_eos_deriv(1,1) = detRow(nr_eos_k)/abs(nr_eos_eint)*nr_eos_scale(1)
        elseif (nr_eos_mode .eq. MODE_PRES_ENTR) then
            newt_eos(1) = (ptotRow(nr_eos_k) - nr_eos_pres)/abs(nr_eos_pres)
            newt_eos(2) = (stotRow(nr_eos_k) - nr_eos_entr)/abs(nr_eos_entr)
            nr_eos_deriv(1,1) = dpdRow(nr_eos_k)/abs(nr_eos_pres)*nr_eos_scale(1)
            nr_eos_deriv(1,2) = dptRow(nr_eos_k)/abs(nr_eos_pres)*nr_eos_scale(2)
            nr_eos_deriv(2,1) = dsdRow(nr_eos_k)/abs(nr_eos_entr)*nr_eos_scale(1)
            nr_eos_deriv(2,2) = dstRow(nr_eos_k)/abs(nr_eos_entr)*nr_eos_scale(2)
        elseif (nr_eos_mode .eq. MODE_PRES_TEMP) then
            newt_eos(1) = (ptotRow(nr_eos_k) - nr_eos_pres)/abs(nr_eos_pres)
            nr_eos_deriv(1,1) = dpdRow(nr_eos_k)/abs(nr_eos_pres)*nr_eos_scale(1)
        else
            call Driver_abortFlash('ERROR: Incorrect EOS mode specified.')
        endif
#ifdef DEBUG_EOS
        !if (nr_eos_deriv .eq. 0.0) call Driver_abortFlash('Error: Eos Derivative is zero!')
        if (newt_eos(1) .ne. newt_eos(1)) then
            print *, 'input', tlim
            print *, 'var info', newt_eos(1), denRow(nr_eos_k), ptotRow(nr_eos_k), etotRow(nr_eos_k), nr_eos_pres, nr_eos_eint
            print *, 'deriv info', nr_eos_deriv(1,1), nr_eos_scale(1), nr_eos_mode
            call Driver_abortFlash('Error: Eos is NaN!')
        endif
        if (nr_eos_deriv(1,1) .ne. nr_eos_deriv(1,1)) then
            print *, 'input', tlim
            print *, 'var info', newt_eos(1), ptotRow(nr_eos_k), etotRow(nr_eos_k), nr_eos_pres, nr_eos_eint
            print *, 'deriv info', nr_eos_deriv, nr_eos_scale(1), nr_eos_mode
            call Driver_abortFlash('Error: Eos Derivative is NaN!')
        endif
        if (minval(abs(nr_eos_deriv(1:size(x),1:size(x)))) .eq. 0.0) then
            !newt_eos(1) = 0.0
            call Driver_abortFlash('nr_eos_deriv zero!')
        else
#endif
            if (nr_eos_limited) then
                !if (nr_eos_mode .eq. MODE_DENS_EI) then
                !    print *, nr_eos_scale, newt_eos, nr_eos_deriv, x
                !    call Driver_abortFlash('limited')
                !endif
                nr_eos_lim_cnt = nr_eos_lim_cnt + 1
                if (nr_eos_lim_cnt .gt. 10 .or. nr_eos_mode .eq. MODE_PRES_ENTR) then !Hit limit too many times, exit
                    newt_eos(1) = 0.0
                    return
                endif
                if (nr_eos_mode .ne. MODE_PRES_TEMP) then
                    !newt_eos(1) = 0
                    sclim = tlim/nr_eos_scale(1)
                elseif (nr_eos_mode .eq. MODE_PRES_TEMP) then
                    !newt_eos(1) = 0
                    sclim = dlim/nr_eos_scale(1)
                endif
                if (abs(x(1)) .lt. sclim) then
                    newt_eos(1) = newt_eos(1)*log(exp(1.0) + 1.0 - x(1)/sclim)
                    nr_eos_deriv = newt_eos(1)/sclim/(x(1)/sclim - exp(1.0) - 1.0)
                else
                    newt_eos(1) = newt_eos(1)*log(exp(1.0) + x(1)/sclim - 1.0)
                    nr_eos_deriv = newt_eos(1)/sclim/(x(1)/sclim + exp(1.0) - 1.0)
                endif
#ifdef DEBUG_EOS
            endif
#endif
        endif
    END FUNCTION newt_eos

END MODULE eos_newt_functions
