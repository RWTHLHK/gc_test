//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html
#include "gc_testTestApp.h"
#include "gc_testApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "MooseSyntax.h"

InputParameters
gc_testTestApp::validParams()
{
  InputParameters params = gc_testApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  return params;
}

gc_testTestApp::gc_testTestApp(InputParameters parameters) : MooseApp(parameters)
{
  gc_testTestApp::registerAll(
      _factory, _action_factory, _syntax, getParam<bool>("allow_test_objects"));
}

gc_testTestApp::~gc_testTestApp() {}

void
gc_testTestApp::registerAll(Factory & f, ActionFactory & af, Syntax & s, bool use_test_objs)
{
  gc_testApp::registerAll(f, af, s);
  if (use_test_objs)
  {
    Registry::registerObjectsTo(f, {"gc_testTestApp"});
    Registry::registerActionsTo(af, {"gc_testTestApp"});
  }
}

void
gc_testTestApp::registerApps()
{
  registerApp(gc_testApp);
  registerApp(gc_testTestApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
// External entry point for dynamic application loading
extern "C" void
gc_testTestApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  gc_testTestApp::registerAll(f, af, s);
}
extern "C" void
gc_testTestApp__registerApps()
{
  gc_testTestApp::registerApps();
}
