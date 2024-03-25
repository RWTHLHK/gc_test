#include "gc_testApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
gc_testApp::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

gc_testApp::gc_testApp(InputParameters parameters) : MooseApp(parameters)
{
  gc_testApp::registerAll(_factory, _action_factory, _syntax);
}

gc_testApp::~gc_testApp() {}

void 
gc_testApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAllObjects<gc_testApp>(f, af, s);
  Registry::registerObjectsTo(f, {"gc_testApp"});
  Registry::registerActionsTo(af, {"gc_testApp"});

  /* register custom execute flags, action syntax, etc. here */
}

void
gc_testApp::registerApps()
{
  registerApp(gc_testApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
gc_testApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  gc_testApp::registerAll(f, af, s);
}
extern "C" void
gc_testApp__registerApps()
{
  gc_testApp::registerApps();
}
