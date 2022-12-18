module treewalk_subs
    integer, parameter :: short = SELECTED_INT_KIND(2)
    integer, parameter :: long = SELECTED_INT_KIND(12)
    integer, parameter :: qp = selected_real_kind(33, 4931)

contains 

      function get_intv(leaf, n, knot) result(pts)
          implicit none
          integer, intent(in) :: n, knot, leaf(:)
          integer :: k1, k2, k3
          real(qp) :: pts(2), t1, t2, t3, t4
           
          ! knot = position of the zero
          ! t1, t2 - upper bound 
          ! t3, t4 - lower bound
          if ((mod(knot,2) == 1).and.(mod(n-knot,2)==0)) then
             t1 = sqrt(3.0_qp)-1.0_qp
             t2 = t1 
             t3 = t1*0.50_qp
             t4 = t3 
          else if  ((mod(knot,2) == 1).and.(mod(n-knot,2)==1)) then
             t1 = sqrt(3.0_qp)-1.0_qp
             t2 = t1*0.50_qp 
             t3 = t2
             t4 = t1 
          else if  ((mod(knot,2) == 0).and.(mod(n-knot,2)==0)) then
             t2 = sqrt(3.0_qp)-1.0_qp
             t1 = t2*0.50_qp 
             t3 = t2
             t4 = t1 
          else if  ((mod(knot,2) == 0).and.(mod(n-knot,2)==1)) then
             t1 = 0.50_qp*(sqrt(3.0_qp)-1.0_qp)
             t2 = t1 
             t3 = sqrt(3.0_qp)-1.0_qp
             t4 = t3 
          end if  

          if ((n>1).and.(knot>1)) then 
!          x1 = t3; 
            do k1 = 1, knot-1
                  t3 = 1.0_qp/(real(leaf(k1),qp) + t3) 
                  t1 = 1.0_qp/(real(leaf(k1),qp) + t1) 
            enddo
            do k1 = n, knot + 1, - 1
                 t4 = 1.0_qp/(real(leaf(k1),qp) + t4);
                t2 = 1.0_qp/(real(leaf(k1),qp) + t2);
            enddo  
          elseif ((n > 1).and.(knot == 1)) then
               do k1 = n, knot + 1, - 1
                 t4 = 1.0_qp/(real(leaf(k1),qp) + t4);
                 t2 = 1.0_qp/(real(leaf(k1),qp) + t2);
               enddo
          endif
          
          pts(1) = t4 + t3 + leaf(knot)
          pts(2) = t2 + t1 + leaf(knot)  

      end function get_intv

      recursive subroutine fruits(leaf, n, knot, t0, outfile, leftb, maxit, ct, lm)
          implicit none
          integer, intent(in) :: maxit, n, knot
          real(qp), intent(in) :: t0
          integer, intent(inout) :: leaf(:), ct, lm
          real(qp), intent(inout) :: leftb

          real(qp) :: pts(2)
          integer :: j
          integer :: leaf_copy(maxit)
         
          character(100), intent(in) :: outfile
          character(25) :: myprint, dimchar

          leaf_copy = 0

          if (n > 0) then
            pts = get_intv(leaf,n,knot)

            if (pts(1) > t0) then  
                open(7, file = trim(outfile), status = 'old', position='append')
                write (myprint, '(A1I2A7)') trim("("), n+1, trim("(1xI3))")
                write (7, FMT=myprint) n, (leaf(j), j=1, n )
                write (7, FMT=myprint) n, (leaf(j), j=n,1,-1 )
                !print *, n, (leaf(j), j=1, n )
                !print *, n, (leaf(j), j=n,1,-1 )
                ct = ct + 2
                if (lm < n) then 
                    lm = n
                endif
                close(7)
            elseif (pts(2) < t0) then 
                if (pts(2) > leftb) then 
                      leftb = pts(2)
                endif
            elseif (n < maxit) then
                if (mod(n,2) == 0) then
                  leaf(n+1) = 1
                  leaf_copy = leaf
                  call fruits(leaf,n+1,knot,t0,outfile,leftb,maxit,ct,lm)
                  leaf_copy(n+1) = 2 
                  call fruits(leaf_copy,n+1,knot,t0,outfile,leftb,maxit,ct,lm)
                else  
                  do j = n+1, 2, -1
                    leaf(j) = leaf(j-1)
                  enddo
                    leaf(1) = 1
                    leaf_copy = leaf
                    call fruits(leaf,n+1,knot+1,t0,outfile,leftb,maxit,ct,lm)
                    leaf_copy(1) = 2
                    call fruits(leaf_copy,n+1,knot+1,t0,outfile,leftb,maxit,ct,lm)
                endif
            elseif (pts(2) > leftb) then 
                   leftb = pts(2)
                !   print *, 'Last step, the leaf would be divided: n=', n
                !endif
            endif
          else 
             leaf(1) = 1 
             call fruits(leaf,1,1,t0,outfile,leftb,maxit,ct,lm)
             leaf_copy(1) = 2    
             call fruits(leaf_copy,1,1,t0,outfile,leftb,maxit,ct,lm)
         endif
            
      end subroutine fruits  

endmodule treewalk_subs

program treewalk
    use treewalk_subs
      implicit none
      integer, allocatable :: leaf(:)
      integer :: n, knot, maxit, num, ct, lm
      real(qp) :: leftb, t0 
      character(100) :: outfile, dimchar

     ! maxit = 7
      n = 0 
      ct = 0
      lm = 0
     ! t0 = 3.10_qp 
      leftb = 0.0_qp 

      if(command_argument_count()<3) then
        write(*,*) 'The value t_0, the search level and file number are required.' 
        stop
      else   
        call get_command_argument(1,dimchar) 
        read(dimchar,*) t0
        call get_command_argument(2,dimchar) 
        read(dimchar,*) maxit
        call get_command_argument(3,dimchar) 
        read(dimchar,*) num
      endif 

      allocate(leaf(maxit))
      write (outfile, '(A6I4A4)') trim("fruits"), num, trim(".dat")
      open(12, file = trim(outfile))
      write(12,*) '2'
      write(12,*) '1     2'
      close(12)
     ! print *, 't0:', t0, 'looking for sequences up to length', maxit

      call fruits(leaf, n, knot, t0, outfile, leftb, maxit, ct, lm)
      print *, 't0:', t0, 'lb:', leftb
      open(12, file = trim(outfile), status = 'old', position='append')
            write (12, '(I3)')  lm-1
            write (12, '(I3)')  ct
            write (12,*)  '4'  
      close(12)

      deallocate(leaf)


end program treewalk
           
