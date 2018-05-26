       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Tue Dec 19 13:08:05 2017'
       b_stamp_str = 'Fri Dec 22 10:03:56 2017'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& login2.stampede.tacc.utexas.edu&
& 2.6.32-431.17.1.el6.x86_64&
& #1 SMP Wed May 7 23:32:49 UTC 2014&
& x86_64'
       return
       end subroutine
      
