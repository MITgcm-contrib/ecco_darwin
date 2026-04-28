      subroutine create_directory(path)
          use iso_c_binding
          implicit none
          character*(*) path

          interface
              integer function mkdir(cpath, mode) bind(C, name="mkdir")
                  use iso_c_binding
                  character(kind=c_char), dimension(*) :: cpath
                  integer(c_int16_t), value :: mode
              end function mkdir
          end interface

          integer :: status
          character(len=100) :: cpath
          integer :: i

          ! Convert Fortran string to C string
          do i = 1, len_trim(path)
              cpath(i:i) = path(i:i)
          end do
          cpath(len_trim(path)+1:len_trim(path)+1) = char(0)

          status = mkdir(cpath, int(o'755',c_int16_t))

      end subroutine create_directory

