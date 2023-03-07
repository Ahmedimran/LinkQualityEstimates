/*
 * MAPPLE-online-app-helper.h
 *
 *  Created on: Jul 18, 2012
 *      Author: imran
 */

#ifndef MAPPLE_ONLINE_APP_HELPER_H_
#define MAPPLE_ONLINE_APP_HELPER_H_


#include "ns3/object-factory.h"
#include "ns3/ipv4-address.h"
#include "ns3/node-container.h"
#include "ns3/application-container.h"
#include "ns3/ipv4-interface-container.h"

namespace ns3 {

/**
 * \brief A helper to make it easier to instantiate an ns3::PacketSinkApplication
 * on a set of nodes.
 */
class MappleOnlineApplicationkHelper
{
public:
  /**
   * Create a PacketSinkHelper to make it easier to work with PacketSinkApplications
   *
   * \param protocol the name of the protocol to use to receive traffic
   *        This string identifies the socket factory type used to create
   *        sockets for the applications.  A typical value would be
   *        ns3::TcpSocketFactory.
   * \param address the address of the sink,
   *
   */
	MappleOnlineApplicationkHelper (std::string protocol, uint16_t portBase, Ipv4InterfaceContainer inter, NodeContainer nc, std::string prefix,double simtime);

  /**
   * Helper function used to set the underlying application attributes.
   *
   * \param name the name of the application attribute to set
   * \param value the value of the application attribute to set
   */
  void SetAttribute (std::string name, const AttributeValue &value);

  /**
   * Install an ns3::PacketSinkApplication on each node of the input container
   * configured with all the attributes set with SetAttribute.
   *
   * \param c NodeContainer of the set of nodes on which a PacketSinkApplication
   * will be installed.
   */
  ApplicationContainer Install (NodeContainer c) const;

  /**
   * Install an ns3::PacketSinkApplication on each node of the input container
   * configured with all the attributes set with SetAttribute.
   *
   * \param node The node on which a PacketSinkApplication will be installed.
   */
  ApplicationContainer Install (Ptr<Node> node) const;

  /**
   * Install an ns3::PacketSinkApplication on each node of the input container
   * configured with all the attributes set with SetAttribute.
   *
   * \param nodeName The name of the node on which a PacketSinkApplication will be installed.
   */
  ApplicationContainer Install (std::string nodeName) const;

private:
  /**
   * \internal
   */
  Ptr<Application> InstallPriv (Ptr<Node> node) const;
  ObjectFactory m_factory;

  Ipv4InterfaceContainer iface;
  NodeContainer nodes;
  std::string prefix;
  double simulationtime;
  //std::string model;
};

} // namespace ns3


#endif /* MAPPLE_ONLINE_APP_HELPER_H_ */
