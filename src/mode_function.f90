integer function stat_mode(a, others, otherslen, ok)
    integer, dimension(:), intent(in) :: a
    logical, optional, intent(out)    :: ok
    integer, dimension(size(a,1)), optional, intent(out) :: others
    integer, optional, intent(out)    :: otherslen

    ! ta is a copy of a, we sort ta modifying it, freq
    ! holds the frequencies and idx the index (for ta) so that
    ! the value appearing freq(i)-time is ta(idx(i))
    integer, dimension(size(a, 1)) :: ta, freq, idx
    integer                        :: rs, i, tm, ml, tf

    if ( present(ok) ) ok = .false.

    select case ( size(a, 1) )
    case (0)  ! no mode... ok is false
       return
    case (1)
       if ( present(ok) ) ok = .true.
       stat_mode = a(1)
       return
    case default
       if ( present(ok) ) ok = .true.
       ta = a         ! copy the array
       call sort(ta)  ! sort it in place (cfr. sort algos on RC)
       freq = 1
       idx = 0
       rs = 1         ! rs will be the number of different values
       
       do i = 2, size(ta, 1)
          if ( ta(i-1) == ta(i) ) then
             freq(rs) = freq(rs) + 1
          else
             idx(rs) = i-1
             rs = rs + 1
          end if
       end do
       idx(rs) = i-1
         
       ml = maxloc(freq(1:rs), 1)  ! index of the max value of freq
       tf = freq(ml)               ! the max frequency
       tm = ta(idx(ml))            ! the value with that freq

       ! if we want all the possible modes, we provide others
       if ( present(others) ) then
          i = 1
          others(1) = tm
          do
             freq(ml) = 0
             ml = maxloc(freq(1:rs), 1)
             if ( tf == freq(ml) ) then ! the same freq
                i = i + 1               ! as the max one
                others(i) = ta(idx(ml))
             else
                exit
             end if
          end do
                
          if ( present(otherslen) ) then
             otherslen = i
          end if

       end if
       stat_mode = tm
    end select

  end function stat_mode
