/*
 * MAPPLE-online-app-helper.cc
 *
 *  Created on: Jul 18, 2012
 *      Author: imran
 */


#include "MAPPLE-online-app-helper.h"
#include "ns3/MAPPLE-online-app.h"
#include "ns3/string.h"
#include "ns3/inet-socket-address.h"
#include "ns3/names.h"
#include "ns3/uinteger.h"

namespace ns3 {

MappleOnlineApplicationkHelper::MappleOnlineApplicationkHelper (std::string protocol, uint16_t portBase, Ipv4InterfaceContainer inter, NodeContainer nc, std::string prefix,double simtime)
{
  m_factory.SetTypeId ("ns3::MappleOnlineApplication");
  m_factory.Set ("Protocol", StringValue (protocol));
  m_factory.Set("port", UintegerValue(portBase));

  this -> iface = inter;
  this -> nodes = nc;
  this -> prefix = prefix;
  this -> simulationtime = simtime;
  //this -> model = model;
}

void
MappleOnlineApplicationkHelper::SetAttribute (std::string name, const AttributeValue &value)
{
  m_factory.Set (name, value);
}

ApplicationContainer
MappleOnlineApplicationkHelper::Install (Ptr<Node> node) const
{
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
MappleOnlineApplicationkHelper::Install (std::string nodeName) const
{
  Ptr<Node> node = Names::Find<Node> (nodeName);
  return ApplicationContainer (InstallPriv (node));
}

ApplicationContainer
MappleOnlineApplicationkHelper::Install (NodeContainer c) const
{
  ApplicationContainer apps;
  for (NodeContainer::Iterator i = c.Begin (); i != c.End (); ++i)
    {
      apps.Add (InstallPriv (*i));
    }

  return apps;
}

Ptr<Application>
MappleOnlineApplicationkHelper::InstallPriv (Ptr<Node> node) const
{
	//Refference to our application using Base class pointer refference which is Application///
  Ptr<Application> app = m_factory.Create<Application> ();
  Ptr<MappleOnlineApplication> mappleApp = DynamicCast<MappleOnlineApplication>(app);
  mappleApp->SetInterfaceContainer(this->iface);
  mappleApp->SetSimulationtime(this->simulationtime);
  mappleApp->SetNodeContainer(this->nodes);
  if ( (this->prefix) != "" )
	  mappleApp->SetPrefix(this->prefix);
//  if ( (this->model) != "" )
//	  mappleApp->SetModel(this->model);
  node->AddApplication (app);

  return app;
}

} // namespace ns3




