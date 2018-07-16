MODULE BRUSSELATOR_IO
    use mpi
    use decomp_2d
    use adios2
    public 

    type(adios2_adios)      :: adios2_handle
    type(adios2_io)         :: io_obj
    type(adios2_engine)     :: engine
    type(adios2_variable)   :: var_plotnum, var_field
    logical                 :: adios2_initialized

    CONTAINS

    SUBROUTINE savedata(Nx,Ny,Nz,plotnum,name,field,u,v,decomp)
        !--------------------------------------------------------------------
        !
        !
        ! PURPOSE
        !
        ! This subroutine saves a three dimensional real array in binary 
        ! format
        !
        ! INPUT
        !
        ! .. Scalars ..
        !  Nx				= number of modes in x - power of 2 for FFT
        !  Ny				= number of modes in y - power of 2 for FFT
        !  Nz				= number of modes in z - power of 2 for FFT
        !  plotnum			= number of plot to be made
        ! .. Arrays ..
        !  field 			= real data to be saved
        !  name_config		= root of filename to save to 
        !
        ! .. Output	..	
        ! plotnum			= number of plot to be saved
        ! .. Special Structures ..
        !  decomp			= contains information on domain decomposition
        !					see http://www.2decomp.org/ for more information
        ! LOCAL VARIABLES
        !
        ! .. Scalars ..
        !  i				= loop counter in x direction
        !  j				= loop counter in y direction
        !  k				= loop counter in z direction
        !  count			= counter
        !  iol				= size of file
        ! .. Arrays ..
        ! 	number_file		= array to hold the number of the plot
        !
        ! REFERENCES
        !
        ! ACKNOWLEDGEMENTS
        !
        ! ACCURACY
        !		
        ! ERROR INDICATORS AND WARNINGS
        !
        ! FURTHER COMMENTS
        !--------------------------------------------------------------------
        ! External routines required
        ! 
        ! External libraries required
        ! 2DECOMP&FFT	 -- Domain decomposition and Fast Fourier Library
        !			(http://www.2decomp.org/index.html)
        ! MPI library
        USE decomp_2d
        USE decomp_2d_fft
        USE decomp_2d_io
        IMPLICIT NONE					 
        INCLUDE 'mpif.h'
        ! Declare variables
        INTEGER(KIND=4), INTENT(IN)						:: Nx,Ny,Nz
        INTEGER(KIND=4), INTENT(IN)						:: plotnum
        TYPE(DECOMP_INFO), INTENT(IN)					::  decomp
        REAL(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
            decomp%xst(2):decomp%xen(2),&
            decomp%xst(3):decomp%xen(3)), &
            INTENT(INOUT) :: field
        COMPLEX(KIND=8), DIMENSION(decomp%xst(1):decomp%xen(1),&
            decomp%xst(2):decomp%xen(2),&
            decomp%xst(3):decomp%xen(3)), &
            INTENT(IN) :: u,v
        CHARACTER*100, INTENT(IN)	     				:: name
        CHARACTER*100               					:: name_config
        INTEGER(kind=4)									:: i,j,k,iol,count,ind
        CHARACTER*100									:: number_file
    
#ifdef ADIOS2
        call savedata_adios2(Nx,Ny,Nz,plotnum,name,field,u,v,decomp)
#else
    
        ! create character array with full filename
        ! write out using 2DECOMP&FFT MPI-IO routines
        ind=index(name,' ') -1
        name_config=name(1:ind)//'u'
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        field(i,j,k)=REAL(u(i,j,k))
        END DO; END DO; END DO
        ind = index(name_config,' ') - 1
        WRITE(number_file,'(i0)') plotnum
        number_file = name_config(1:ind)//number_file
        ind = index(number_file,' ') - 1
        number_file = number_file(1:ind)//'.datbin'	
        CALL decomp_2d_write_one(1,field,number_file, decomp)
    
        ind=index(name,' ') -1
        name_config=name(1:ind)//'v'
        DO k=decomp%xst(3),decomp%xen(3); DO j=decomp%xst(2),decomp%xen(2); DO i=decomp%xst(1),decomp%xen(1)
        field(i,j,k)=REAL(v(i,j,k))
        END DO; END DO; END DO
        ind = index(name_config,' ') - 1
        WRITE(number_file,'(i0)') plotnum
        number_file = name_config(1:ind)//number_file
        ind = index(number_file,' ') - 1
        number_file = number_file(1:ind)//'.datbin'	
        CALL decomp_2d_write_one(1,field,number_file, decomp)
    
#endif
    END SUBROUTINE savedata
    
    
    !----------------------------------------------------------------------------!
    SUBROUTINE io_init (decomp, ierr)
    !----------------------------------------------------------------------------!
        implicit none
        
        type(decomp_info), intent(in)   :: decomp
        integer, intent(out)            :: ierr
        integer*8, dimension(3)         :: sizes, subsizes, starts

        ierr = 0
#ifdef ADIOS2
        ! Code for getting sizes, subsizes, and starts copied from 2decomp_fft
        ! determine subarray parameters
        sizes(1) = decomp%xsz(1)
        sizes(2) = decomp%ysz(2)
        sizes(3) = decomp%zsz(3)
        
        subsizes(1) = decomp%xsz(1)
        subsizes(2) = decomp%xsz(2)
        subsizes(3) = decomp%xsz(3)
        starts(1) = decomp%xst(1)-1  ! 0-based index
        starts(2) = decomp%xst(2)-1
        starts(3) = decomp%xst(3)-1
        
        ! Init adios2
        call adios2_init_config (adios2_handle, "adios2_config.xml", mpi_comm_world, &
            adios2_debug_mode_off, ierr)
    
        ! Init IO object
        call adios2_declare_io (io_obj, adios2_handle, "savedata", ierr)
    
        ! Define variables
        call adios2_define_variable (var_plotnum, io_obj, "plotnum", adios2_type_integer4, ierr)
        call adios2_define_variable (var_field, io_obj, "field", adios2_type_dp, 3, &
            sizes, starts, subsizes, .true., ierr)
        
        ! Open file
        call adios2_open (engine, io_obj, "./data/brusselator", adios2_mode_write, &
            mpi_comm_world, ierr)
#endif
    end subroutine io_init
    !----------------------------------------------------------------------------!


    !----------------------------------------------------------------------------!
    SUBROUTINE savedata_adios2(Nx,Ny,Nz,plotnum,name,field,u,v,decomp)
    !----------------------------------------------------------------------------!
        use decomp_2d
        use decomp_2d_fft
        implicit none					 
        include 'mpif.h'
        integer(kind=4), intent(in)						                :: Nx,Ny,Nz
        integer(kind=4), intent(in)						                :: plotnum
        type(decomp_info), intent(in)					                ::  decomp
        real(kind=8), dimension(decomp%xst(1):decomp%xen(1),&
                                decomp%xst(2):decomp%xen(2),&
                                decomp%xst(3):decomp%xen(3)), &
                      intent(inout)                                     :: field
        complex(kind=8), dimension(decomp%xst(1):decomp%xen(1),&
                                   decomp%xst(2):decomp%xen(2),&
                                   decomp%xst(3):decomp%xen(3)), &
                         intent(in)                                     :: u,v
        
        character*100, intent(in)	     				                :: name
        character*100               					                :: name_config
        integer(KIND=4)									                :: i,j,k,iol,count,ind
        character*100									                :: number_file
        integer                                                         :: myrank, ierr
    
        call mpi_comm_rank (mpi_comm_world, myrank, ierr)

        call adios2_begin_step  (engine, ierr)
        if (myrank .eq. 0) then
            call adios2_put     (engine, var_plotnum, plotnum, ierr)
        endif
        call adios2_put         (engine, var_field, field, ierr)
        call adios2_end_step    (engine, ierr)
    
    END SUBROUTINE savedata_adios2
    
    
    !----------------------------------------------------------------------------!
    SUBROUTINE IO_FINALIZE(ierr)
    !----------------------------------------------------------------------------!
        implicit none
    
        integer, intent(out) :: ierr
    
        call adios2_close(engine, ierr)
    
    END SUBROUTINE IO_FINALIZE

END MODULE BRUSSELATOR_IO

