subroutine RadTrans_getDbgContext(context)
  use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
  use RadTrans_data, ONLY: rt_dbgContext
  implicit none
  type(RadTrans_dbgContext_t),intent(OUT) :: context
  
  context = rt_dbgContext

end subroutine RadTrans_getDbgContext

subroutine RadTrans_getDbgContextPtr(context)
  use RadTrans_interfaceTypeDecl, ONLY: RadTrans_dbgContext_t
  use RadTrans_data, ONLY: rt_dbgContext
  implicit none
  type(RadTrans_dbgContext_t),pointer :: context
  
  context => rt_dbgContext

end subroutine RadTrans_getDbgContextPtr
