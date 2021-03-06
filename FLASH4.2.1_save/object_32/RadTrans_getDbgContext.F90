subroutine RadTrans_getDbgContext(context)
  use RadTrans_interface, ONLY: RadTrans_dbgContext_t
  implicit none
  type(RadTrans_dbgContext_t),intent(OUT) :: context
  
  context%step = 0
  context%group = 0

end subroutine RadTrans_getDbgContext

subroutine RadTrans_getDbgContextPtr(context)
  use RadTrans_interface, ONLY: RadTrans_dbgContext_t
  implicit none
  type(RadTrans_dbgContext_t),pointer :: context
  
  nullify(context)

end subroutine RadTrans_getDbgContextPtr
