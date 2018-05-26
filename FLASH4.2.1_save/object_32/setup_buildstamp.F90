       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Mon Oct 26 12:44:59 2015'
       b_stamp_str = 'Mon Oct 26 14:30:54 2015'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& login1.stampede.tacc.utexas.edu&
& 2.6.32-431.17.1.el6.x86_64&
& #1 SMP Wed May 7 23:32:49 UTC 2014&
& x86_64'
       return
       end subroutine
      
