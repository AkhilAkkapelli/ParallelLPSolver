module projective_scaling_mod
use mpi
use mkl_vsl
implicit none

integer, parameter :: m  = 2**8, n = 2**9
integer, parameter :: mb = 2**5, nb = 2**6
integer :: nprow = 3, npcol = 4

integer, parameter :: iterlimit = 100

character(len=8), dimension(4), parameter :: colors = [ char(27)//"[31m", char(27)//"[32m", char(27)//"[33m", char(27)//"[34m" ]

integer, external :: numroc
integer, external :: indxl2g, indxg2l

contains

  subroutine init_mpi_blacs(myrow, mycol, rank, nprocs, ictxt)

    integer, intent(out) :: myrow, mycol
    integer, intent(out) :: rank, nprocs
    integer, intent(out) :: ictxt

    integer :: info

    call MPI_Init(info)

    call blacs_pinfo(rank, nprocs)
    call blacs_get(-1, 0, ictxt)

    if (nprocs /= nprow * npcol) then
    if (rank == 0) print *, "ERROR: This program requires exactly 4 MPI ranks"
      call MPI_Abort(MPI_COMM_WORLD, 1, info)
    end if

    call blacs_gridinit(ictxt,   'R', nprow, npcol)
    call blacs_gridinfo(ictxt, nprow, npcol, myrow, mycol)

    print('(A,I0,A,I0,A,I0,A,I0,A)'), 'BLACS rank:', rank, '/', nprocs, ' BLACS Grid position: (', myrow, ',', mycol, ')'

  end subroutine init_mpi_blacs

  subroutine init_lp(A, x0, b, c, mloc, nloc, descA, descx0, descb, descc, myrow, mycol, ictxt)

    real(8), allocatable, intent(out) :: A(:,:), x0(:,:), b(:,:), c(:,:)
    integer, intent(out) :: mloc, nloc
    integer, intent(out) :: descA(9), descx0(9), descb(9), descc(9)

    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: ictxt

    integer :: info

    mloc = numroc(m, mb, myrow, 0, nprow)
    nloc = numroc(n, nb, mycol, 0, npcol)

    allocate(A (mloc, nloc))
    allocate(x0(   1, nloc))
    allocate(b (mloc,    1))
    allocate(c (   1, nloc))

    call descinit( descA,    m,   n,   mb, nb, 0, 0, ictxt, mloc, info)
    call descinit(descx0,    1,   n,    1, nb, 0, 0, ictxt,    1, info)
    call descinit( descb,    m,   1,   mb,  1, 0, 0, ictxt, mloc, info)
    call descinit( descc,    1,   n,    1, nb, 0 ,0 ,ictxt,    1, info)

    print '(A,I0,A,I0, A,I0,A,I0,A, A,I0,A,I0,A, A,I0,A, A,I0,A, A,I0,A)', 'Rank ', rank, '/', nprocs, &
            ' BLACS(', myrow, ',', mycol, ')', ' A(', mloc, ',', nloc, ')', ' x0(', nloc, ')', ' b(', mloc, ')', ' c(', nloc, ')'

    call print_desc(descA,  'A')
    call print_desc(descx0, 'x0')
    call print_desc(descb,  'b')
    call print_desc(descc,  'c')

  end subroutine init_lp

  subroutine print_desc(desc, name, rank)
    integer, intent(in) :: desc(9)
    character(len=*), intent(in) :: name
    integer, intent(in) :: rank

    if (rank == 0) then
      print ('(A)'), '==========================================='
      print ('(A)'), 'Descriptor for '//trim(name)//':'
      print ('(A)'), '-------------------------------------------'
      print ('(A,I0,A,I0)'), 'Global size        : ', desc(3), ' x ', desc(4)
      print ('(A,I0,A,I0)'), 'Blocking factors   : ', desc(5), ' x ', desc(6)
      print ('(A,I0,A,I0)'), 'Source process     : ', desc(7), ' , ', desc(8)
      print ('(A,I0)'), 'Local leading dim  : ', desc(9)
      print ('(A)'), '-------------------------------------------'
    end if

  end subroutine print_desc

  subroutine gen_data(A, x0, b, c, descA, descx0, descb, mloc, nloc, myrow, mycol, rank, ictxt)

    real(8), intent(inout) :: A(:,:), x0(:,:), b(:,:), c(:,:)

    integer, intent(in) :: descA(9), descx0(9), descb(9)
    
    integer, intent(in) :: mloc, nloc
    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: rank
    integer, intent(in) :: ictxt

    integer :: vsl_seed
    integer :: vsl_status
    type(VSL_STREAM_STATE) :: stream

    vsl_seed   = 12345 + rank * 10
    vsl_status = vslnewstream(stream, VSL_BRNG_MT19937, vsl_seed)

    if (mloc*nloc > 0) then
      vsl_status = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, mloc*nloc, A, 0.0d0, 1.0d0 )
    end if

    if (myrow == 0 .and. nloc > 0) then
      vsl_status = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, nloc, x0, 0.0d0, 1.0d0 )
      x0 = abs(x0)

      vsl_status = vdrnggaussian( VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, nloc, c, 0.0d0, 1.0d0 )
    end if

    vsl_status = vsldeletestream(stream)

    if (myrow == 0) then
      call DGEBS2D(ictxt, 'C', ' ', nloc, 1, x0, nloc)
      call DGEBS2D(ictxt, 'C', ' ', nloc, 1, c,  nloc)
    else
      call DGEBR2D(ictxt, 'C', ' ', nloc, 1, x0, nloc, 0, mycol)
      call DGEBR2D(ictxt, 'C', ' ', nloc, 1, c,  nloc, 0, mycol)
    end if

    call pdgemv( 'N', m, n, 1.0d0, A, 1, 1, descA, x0, 1, 1, descx0, 1, 0.0d0, b, 1, 1, descb, 1 )

    if (mycol == 0) then
      call DGEBS2D(ictxt, 'R', ' ', mloc, 1, b, mloc)
    else
      call DGEBR2D(ictxt, 'R', ' ', mloc, 1, b, mloc, myrow, 0)
    end if

    do p = 0, nprocs - 1
      call MPI_Barrier(MPI_COMM_WORLD, info)

      if (rank == p) then

        print '(A,I0,A,I0,A)', 'BLACS(', myrow, ',', mycol, ')'

        print '(A,I0,A,I0,A)', 'A(', mloc, ',', nloc, '):'
        do i = 1, mloc
          do j = 1, nloc
            write(*,'(A,F10.7,1X)', advance='no') colors(mod( ((i-1)/2 + (j-1)/2), size(colors) ) + 1), A(i,j)
          end do
          write(*,'(A)') char(27)//"[0m"
        end do

        print '(A,I0,A)', 'x0(', nloc, '):'
        do j = 1, nloc
          write(*,'(A,F10.7,1X)', advance='no') colors(mod( (j-1)/2, size(colors) ) + 1), x0(1,j)
        end do
        write(*,'(A)') char(27)//"[0m"

        print '(A,I0,A)', 'c(', nloc, '):'
        do j = 1, nloc
          write(*,'(A,F10.7,1X)', advance='no') colors(mod( (j-1)/2, size(colors) ) + 1), c(1,j)
        end do
        write(*,'(A)') char(27)//"[0m"

        print '(A,I0,A)', 'b(', mloc, '):'
        do i = 1, mloc
          write(*,'(A,F10.7,1X)', advance='no') colors(mod( (i-1)/2, size(colors) ) + 1), b(i,1)
        end do
        write(*,'(A)') char(27)//"[0m"

        print *, ""

      end if
    end do

  end subroutine gen_data

  subroutine write_files(A, x0, b, c, mloc, nloc, myrow, mycol, rank)

    real(8), intent(in) :: A(:,:), x0(:,:), b(:,:), c(:,:)
    integer, intent(in) :: mloc, nloc
    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: rank

    integer :: info
    integer :: fh
    character(len=16) :: buf
    integer :: status(MPI_STATUS_SIZE)
    integer :: iloc, jloc, ig, jg
    integer(MPI_OFFSET_KIND) :: offset
    integer :: FW

    FW = 16

    call MPI_File_open(MPI_COMM_WORLD, 'Aeq.dat', MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, info)

    if (rank == 0) then
      write(buf,'(I7,1X,I7)') m, n
      buf(16:16) = new_line('a')
      call MPI_File_write_at(fh, 0_MPI_OFFSET_KIND, buf, len(buf), MPI_CHARACTER, status, info)
    end if

    do jloc = 1, nloc
      jg = indxl2g(jloc, nb, mycol, 0, npcol)
      if (jg < 1 .or. jg > n) cycle

      do iloc = 1, mloc
        ig = indxl2g(iloc, mb, myrow, 0, nprow)
        if (ig < 1 .or. ig > m) cycle

        write(buf,'(F15.8,1X)') A(iloc,jloc)

        offset = 16 + int((ig-1)*(n*16+1) + (jg-1)*16, MPI_OFFSET_KIND)
        call MPI_File_write_at(fh, offset, buf, FW, MPI_CHARACTER, status, info)

        if (jg == n) then
          offset = 16 + int((ig-1)*(n*16+1) + n*16, MPI_OFFSET_KIND)
          call MPI_File_write_at(fh, offset, new_line('a'), 1, MPI_CHARACTER, status, info)
        end if
      end do
    end do

    call MPI_File_close(fh, info)

    call MPI_File_open(MPI_COMM_WORLD, 'x0eq.dat', MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, info)

    if (rank == 0) then
      write(buf,'(I7,1X,I7)') 1, n
      buf(16:16) = new_line('a')
      call MPI_File_write_at( fh, 0_MPI_OFFSET_KIND, buf, len(buf), MPI_CHARACTER, MPI_STATUS_IGNORE, info )
    end if

    if (myrow == 0 .and. nloc > 0) then

      do jloc = 1, nloc
        jg = indxl2g(jloc, nb, mycol, 0, npcol)
        if (jg < 1 .or. jg > n) cycle

        write(buf,'(F15.8,1X)') x0(1,jloc)

        offset = 16 + int((jg-1) * FW, MPI_OFFSET_KIND)

        call MPI_File_write_at( fh, offset, buf, FW, MPI_CHARACTER, MPI_STATUS_IGNORE, info )
      end do

    end if

    call MPI_File_close(fh, info)

    call MPI_File_open(MPI_COMM_WORLD, 'ceq.dat', MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, info)

    if (rank == 0) then
      write(buf,'(I7,1X,I7)') 1, n
      buf(16:16) = new_line('a')
      call MPI_File_write_at( fh, 0_MPI_OFFSET_KIND, buf, len(buf), MPI_CHARACTER, MPI_STATUS_IGNORE, info )
    end if

    if (myrow == 0 .and. nloc > 0) then

      do jloc = 1, nloc
        jg = indxl2g(jloc, nb, mycol, 0, npcol)
        if (jg < 1 .or. jg > n) cycle

        write(buf,'(F15.8,1X)') c(1,jloc)

        offset = 16 + int((jg-1) * FW, MPI_OFFSET_KIND)
        call MPI_File_write_at( fh, offset, buf, FW, MPI_CHARACTER, MPI_STATUS_IGNORE, info )
      end do

    end if

    call MPI_File_close(fh, info)

    call MPI_File_open(MPI_COMM_WORLD, 'beq.dat', MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, info)

    if (rank == 0) then
      write(buf,'(I7,1X,I7)') m, 1
      buf(16:16) = new_line('a')
      call MPI_File_write_at( fh, 0_MPI_OFFSET_KIND, buf, len(buf), MPI_CHARACTER, MPI_STATUS_IGNORE, info )
    end if

    if (mycol == 0 .and. mloc > 0) then

      do iloc = 1, mloc
        ig = indxl2g(iloc, mb, myrow, 0, nprow)
        if (ig < 1 .or. ig > m) cycle

        write(buf,'(F15.8)') b(iloc,1)
        buf(16:16) = new_line('a')

        offset = 16 + int((ig-1) * FW, MPI_OFFSET_KIND)
        call MPI_File_write_at( fh, offset, buf, FW, MPI_CHARACTER, MPI_STATUS_IGNORE, info )
      end do

    end if

    call MPI_File_close(fh, info)

  end subroutine write_files

  subroutine projective_transform(Acan, ccan, ncanloc, descAcan, descccan, A, x0, c, b, mloc, nloc, myrow, mycol, ictxt)

    real(8), allocatable, intent(out) :: Acan(:,:), ccan(:,:)
    integer, intent(out) :: ncanloc
    integer, intent(out) :: descAcan(9), descccan(9)

    real(8), intent(in) :: A(:,:), x0(:,:), c(:,:), b(:,:)
    integer, intent(in) :: mloc, nloc
    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: ictxt

    integer :: info
    integer :: i, j, ig, jg
    real(8) :: ccan_lmax, ccan_gmax

    ncanloc = numroc(n+1, nb, mycol, 0, npcol)

    allocate(Acan(mloc, ncanloc))
    allocate(ccan(   1, ncanloc))

    print '(A,I0,A,I0, A,I0,A,I0,A, A,I0,A,I0,A, A,I0,A)', 'Rank ', rank, '/', nprocs, &
                  ' BLACS(', myrow, ',', mycol, ')', ' Acan(', mloc, ',', ncanloc, ')', ' ccan(', ncanloc, ')'

    call MPI_Barrier(MPI_COMM_WORLD, info)

    call descinit(descAcan, m,   n+1,   mb, nb, 0, 0, ictxt, mloc, info)
    call descinit(descccan, 1,   n+1,    1, nb, 0, 0, ictxt,    1, info)      

    call print_desc(descAcan, "Acan")
    call print_desc(descccan, "ccan")

    do j = 1, nloc
      Acan(:, j) = A(:, j) * x0(1, j)
    end do

    do j = 1, ncanloc
      jg = indxl2g(j, nb, mycol, 0, npcol)
      if (jg <= n) then
        ccan(1, j) = x0(1, j) * c(1, j)
      else 
        Acan(:, j) = -b(:, 1)
        ccan(1, j) = 0.0d0
      end if
    end do  
  
    ccan_lmax = maxval(ccan)
    call MPI_Allreduce(ccan_lmax, ccan_gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)
    ccan = ccan - ccan_gmax

    do p = 0, nprocs - 1
      call MPI_Barrier(MPI_COMM_WORLD, info)

      if (rank == p) then

        print '(A,I0,A,I0,A)', 'BLACS(', myrow, ',', mycol, ')'

        print '(A,I0,A,I0,A)', 'Acan(', mloc, ',', ncanloc, '):'
        do i = 1, mloc
          do j = 1, ncanloc
            write(*,'(A,F10.7,1X)', advance='no') colors(mod( ((i-1)/2 + (j-1)/2), size(colors) ) + 1), Acan(i,j)
          end do
          write(*,'(A)') char(27)//"[0m"
        end do

        print '(A,I0,A)', 'ccan(', ncanloc, '):'
        do j = 1, ncanloc
          write(*,'(A,F10.7,1X)', advance='no') colors(mod( (j-1)/2, size(colors) ) + 1), ccan(1,j)
        end do
        write(*,'(A)') char(27)//"[0m"

        print *, ""

      end if
    end do

  end subroutine projective_transform

  subroutine inverse_projective_transform(xopt, xoptcan, c, descc, x0, nloc, ncanloc, myrow, mycol, rank, nprocs, ictxt)

    real(8), allocatable, intent(out) :: xopt(:,:)
    real(8), intent(in) :: xoptcan(:,:), x0(:,:)
    real(8), intent(in) :: c(:,:)
    integer, intent(in) :: descc(9)
    integer, intent(in) :: nloc, ncanloc
    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: rank
    integer, intent(in) :: nprocs
    integer, intent(in) :: ictxt

    integer :: info
    real(8) :: xoptcan_last
    integer :: x_gidx
    integer :: x_lastcsrc
    integer :: j, p
    real(8) :: dotx
    integer :: idum

    allocate(xopt(1, nloc))

    x_lastcsrc = mod( (n+1)/nb, npcol )

    if (myrow == 0 .and. mycol == x_lastcsrc) xoptcan_last = xoptcan(1, ncanloc)

    if (myrow == 0 .and. mycol == x_lastcsrc) then
        call dgebs2d(ictxt, 'All', ' ', 1, 1, xoptcan_last, 1)
    else
        call dgebr2d(ictxt, 'All', ' ', 1, 1, xoptcan_last, 1, 0, x_lastcsrc)
    endif

    do j = 1, nloc
      xopt(1,j) = x0(1,j)*xoptcan(1,j)/xoptcan_last
    end do

    call pddot(n, dotx, xopt, 1, 1, descc, 1, c, 1, 1, descc, 1)

    do p = 0, nprocs - 1
      call MPI_Barrier(MPI_COMM_WORLD, info)

      if (rank == p) then

        print '(A,I0,A,I0,A)', 'BLACS(', myrow, ',', mycol, ')'

        print '(A,I0,A)', 'xopt(', nloc, '):'
        do j = 1, nloc
          write(*,'(A,F10.7,1X)', advance='no') colors(mod( (j-1)/2, size(colors) ) + 1), xopt(1,j)
        end do
        write(*,'(A)') char(27)//"[0m"

        print *, ""

      end if
    end do    

    if(myrow == 0 .and. mycol == 0) print*, "Optimum Objective value: ", dotx

  end subroutine inverse_projective_transform

  function projective_scale(Acan, ccan, descccan, ncanloc, mloc, nloc, myrow, mycol, ictxt) result(xoptcan)
 
    real(8), intent(in) :: Acan(:,:), ccan(:,:)
    integer, intent(in) :: descccan(9)
    integer, intent(in) :: ncanloc
    integer, intent(in) :: mloc, nloc
    integer, intent(in) :: myrow, mycol
    integer, intent(in) :: ictxt

    real(8), allocatable :: xoptcan(:,:)
    real(8), allocatable :: xp(:,:), x(:,:)
    integer :: descx(9)
    integer :: info
    integer :: iter
    real(8) :: dotx, dotxp

    iter = 1

    allocate(xp(1, ncanloc))
    allocate( x(1, ncanloc))
    allocate( xoptcan(1, ncanloc))

    call descinit(descx,   1, n+1,  1, nb, 0, 0, ictxt,       1, info)

    xp = 1.0D0/(n+1)

    x = optimize(xp)
    
    if (myrow == 0) then
      call DGEBS2D(ictxt, 'C', ' ', ncanloc, 1, x, ncanloc)
    else
      call DGEBR2D(ictxt, 'C', ' ', ncanloc, 1, x, ncanloc, 0, mycol)
    end if

    call pddot(n+1, dotx, x, 1, 1, descx, 1, ccan, 1, 1, descccan, 1)
    call pddot(n+1, dotxp, xp, 1, 1, descx, 1, ccan, 1, 1, descccan, 1)

    if (myrow == 0 .and. mycol == 0) print*, "iter", iter , dotx

    do while( .NOT. (iter >= iterlimit .OR. dotxp - dotx < 1.0D0-100) )

      iter = iter +1

      xp = x
      x = optimize(xp)
      
      if (myrow == 0) then
        call DGEBS2D(ictxt, 'C', ' ', ncanloc, 1, x, ncanloc)
        call DGEBS2D(ictxt, 'C', ' ', ncanloc, 1, xp, ncanloc)
      else
        call DGEBR2D(ictxt, 'C', ' ', ncanloc, 1, x, ncanloc, 0, mycol)
        call DGEBR2D(ictxt, 'C', ' ', ncanloc, 1, xp, ncanloc, 0, mycol)
      end if

      call pddot(n+1, dotx, x, 1, 1, descx, 1, ccan, 1, 1, descccan, 1)
      call pddot(n+1, dotxp, xp, 1, 1, descx, 1, ccan, 1, 1, descccan, 1)

      if (myrow == 0 .and. mycol == 0) print*, "iter", iter , dotx

    end do

    xoptcan = x

    contains

      function optimize(xp) result(x)

        real(8) ::  x(1, ncanloc)

        real(8), intent(in) ::  xp(:,:)

        integer :: mcanloc
        real(8), allocatable :: B(:,:), r(:,:), BBt(:,:)
        real(8) :: normx, alpha, x_lmax, x_gmax
        integer :: descB(9), descr(9), descBBt(9)
        integer :: i

        mcanloc = numroc(m+1, mb, myrow, 0, nprow)

        allocate(  B(mcanloc, ncanloc))
        allocate(  r(mcanloc,       1))
        allocate(BBt(mcanloc, mcanloc))

        call descinit(   descB, m+1, n+1, mb, nb, 0, 0, ictxt, mcanloc, info)
        call descinit(   descr, m+1,   1, mb,  1, 0, 0, ictxt, mcanloc, info)
        call descinit( descBBt, m+1, m+1, mb, mb, 0, 0, ictxt, mcanloc, info)

        B = 1.0D0
        do i = 1, mloc
          B(i, :) = Acan(i, :) * xp(1, :)
        end do

        x = xp*ccan

        call pdgemv('N', m+1, n+1, 1.0d0, B, 1, 1, descB, x, 1, 1, descx, 1, 0.0d0, r, 1, 1, descr, 1)

        if (mycol == 0) then
          call DGEBS2D(ictxt, 'R', ' ', mcanloc, 1, r, mcanloc)
        else
          call DGEBR2D(ictxt, 'R', ' ', mcanloc, 1, r, mcanloc, myrow, 0)
        end if

        call pdsyrk('L', 'N', m+1, n+1, 1.0d0, B, 1, 1, descB, 0.0d0, BBt, 1, 1, descBBt)

        call pdpotrf('L', m+1, BBt, 1, 1, descBBt, info)
        call pdpotrs('L', m+1, 1, BBt, 1, 1, descBBt, r, 1, 1, descr, info)

        if (mycol == 0) then
          call DGEBS2D(ictxt, 'R', ' ', mcanloc, 1, r, mcanloc)
        else
          call DGEBR2D(ictxt, 'R', ' ', mcanloc, 1, r, mcanloc, myrow, 0)
        end if

        call pdgemv('T', m+1, n+1, 1.0d0, B, 1, 1, descB, r, 1, 1, descr, 1, 0.0d0, x, 1, 1, descx, 1)
        
        x = xp*ccan - x

        call pdnrm2(n+1, normx, x, 1, 1, descx, 1)

        x = x / normx

        x_lmax = maxval(x)
        call MPI_Allreduce(x_lmax, x_gmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, info)

        alpha = 0.1D0 / ( (n+1) * x_gmax )

        x = 1.0D0/(n+1) - alpha*x
        call pddot(n+1, dotx, x, 1, 1, descx, 1, xp, 1, 1, descx, 1)
        
        x = x*xp/dotx
          
      end function optimize

  end function projective_scale

  subroutine cleanup_mpi_blacs(A, x0, b, c, xopt, Acan, ccan, xoptcan, ictxt)

    real(8), allocatable, intent(inout) :: A(:,:), x0(:,:), b(:,:), c(:,:), xopt(:,:), Acan(:,:), ccan(:,:), xoptcan(:,:)
    integer, intent(in) :: ictxt

    deallocate(A, x0, b, c, xopt)
    deallocate(Acan, ccan, xoptcan)

    call blacs_gridexit(ictxt)
    call blacs_exit(0)

  end subroutine cleanup_mpi_blacs

end module projective_scaling_mod

program projective_scaling_algo
    use mpi
    use mkl_vsl
    use projective_scaling_mod
    implicit none

    integer :: rank, nprocs
    integer :: ictxt, myrow, mycol
    integer :: mloc, nloc
    integer :: ncanloc

    real(8), allocatable :: A(:,:), x0(:,:), b(:,:), c(:,:), xopt(:,:)
    real(8), allocatable :: Acan(:,:), ccan(:,:), xoptcan(:,:)

    integer :: descA(9), descx0(9), descb(9), descc(9)
    integer :: descAcan(9), descccan(9)

    integer :: info

    real(8) :: t0, t1
    
    call init_mpi_blacs(myrow, mycol, rank, nprocs, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    call init_lp(A, x0, b, c, mloc, nloc, descA, descx0, descb, descc, myrow, mycol, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    call gen_data(A, x0, b, c, descA, descx0, descb, mloc, nloc, myrow, mycol, rank, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    call write_files(A, x0, b, c, mloc, nloc, myrow, mycol, rank)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    t0 = MPI_Wtime()

    call projective_transform(Acan, ccan, ncanloc, descAcan, descccan, A, x0, c, b, mloc, nloc, myrow, mycol, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    xoptcan = projective_scale(Acan, ccan, descccan, ncanloc, mloc, nloc, myrow, mycol, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)

    call inverse_projective_transform(xopt, xoptcan, c, descc, x0, nloc, ncanloc, myrow, mycol, rank, nprocs, ictxt)

    call MPI_Barrier(MPI_COMM_WORLD, info)
    
    t1 = MPI_Wtime()

    if (rank == 0) print *, "Total time (s): ", t1 - t0

    call cleanup_mpi_blacs(A, x0, b, c, xopt, Acan, ccan, xoptcan, ictxt)

end program projective_scaling_algo
