
#include "ns3/address.h"
#include "ns3/address-utils.h"
#include "ns3/log.h"
#include "ns3/uinteger.h"
#include "ns3/inet-socket-address.h"
#include "ns3/node.h"
#include "ns3/socket.h"
#include "ns3/udp-socket.h"
#include "ns3/simulator.h"
#include "ns3/socket-factory.h"
#include "ns3/packet.h"
#include "ns3/trace-source-accessor.h"
#include "ns3/udp-socket-factory.h"
#include "ns3/wifi-net-device.h"
#include "ns3/regular-wifi-mac.h"
#include "ns3/PhySignalTag.h"
//#include "ns3/MacDelaysTag.h"
//#include "ns3/MacQueueTag.h"
#include "ns3/vector.h"
#include "ns3/ipv4-l3-protocol.h"
#include "ns3/llc-snap-header.h"
#include "ns3/string.h"
#include "ns3/udp-header.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <math.h>
#include <time.h>

#include <inttypes.h>

#include "MAPPLE-online-app.h"

using namespace std;

namespace ns3{

NS_LOG_COMPONENT_DEFINE ("MAPPLEONLINE");
NS_OBJECT_ENSURE_REGISTERED (MappleOnlineApplication);
//
//bool one = false;
//bool two = false;
//bool flag = false;
////Eduardo definitions////

//LWPR_Object lwpr_model(23,1);
FILE * lwpr_file;
FILE* fp_error;
float EWMA_RATE_W = 0.2;
float EMA_RATE_W = 0.013;


UniformVariable uniformvariable(0.0,0.4);


TypeId
MappleOnlineApplication::GetTypeId (void)
{
	static TypeId tid = TypeId ("ns3::MappleOnlineApplication")
			.SetParent<Application> ()
			.AddConstructor<MappleOnlineApplication> ()
			.AddAttribute("port", "Port number on which the ns3 <-> ns3 socket communicates",
			UintegerValue(9),
			MakeUintegerAccessor(&MappleOnlineApplication::portBase),
			MakeUintegerChecker<uint16_t>(1))
			.AddAttribute ("Protocol", "The type id of the protocol to use for the rx socket.",
			TypeIdValue (UdpSocketFactory::GetTypeId ()),
			MakeTypeIdAccessor (&MappleOnlineApplication::m_tid),
			MakeTypeIdChecker ())
			.AddTraceSource ("Rx", "A packet has been received",
			MakeTraceSourceAccessor (&MappleOnlineApplication::m_RecvPacketTrace))
			.AddTraceSource ("Tx", "A new packet is created and is sent",
			MakeTraceSourceAccessor (&MappleOnlineApplication::m_TransPacketTrace))
			;
	return tid;
}
//as we create object of the MappleOnlineApplication it will initiallized the private variables..
MappleOnlineApplication::MappleOnlineApplication ()
{
	NS_LOG_FUNCTION (this);
	m_listningSocket = 0;
	sendingSocket = 0;
	m_summary_counter = 0;
	//	m_my_status = 0;
	m_last_update = 0;
	m_individual_counter = 0;
	last_update = 0;

	oldRate=0;
	oldEstRate=0;
	oldBytesSent=0;

	counter = 0;
	m_Mapp.empty();
	/* Sent packets counter since last summary */
	sent_traffic_counter=0;
	/* Number of samples collected by the node*/
	n_samples = 0;
	ss_counter = 0;
	nooffeatures = 23; //we have exactly 22 features and one 1 target value
	//we statically allocate the memory once, and overwrite the incoming samples again and again..!
	testfeatures = (struct svm_node*) malloc((nooffeatures)*sizeof(struct svm_node));
	//sharedfeatures are those feature which we required for colaborative learning
	//	shared_samples->sharedfeatures = (struct svm_node*) malloc((nooffeatures+2)*sizeof(struct svm_node));

	//recv_sharedsamples = (struct shared_samples_t*) malloc((nooffeatures+2)*sizeof(struct shared_samples_t));

	//oldSamples = (double *) malloc((nooffeatures+2)*sizeof(double));
	//recentSamples = (double *) malloc((nooffeatures+2)*sizeof(double));
	mutex = false;
	check = 1;
	errSum = 0.0;
	errCounter = 1; //intially 1, becuase we post increment

	memset(comparing,0x0,sizeof(double)*24);

	for(uint32_t i=0; i< MAX_NODE_ID;i++)
	{
		memset( &Neighbors[i], 0x0, sizeof(online_neighbor_entry_t));
		memset(&Samples[i], 0x0, sizeof(online_sample_t));
	}


	for(uint32_t i=0; i<5; i++)
		my_status[i] = 0;

	interval =  0;
	starttime = 0;

	hasEMA = 0;
	hasEWMA = 0;
	em_interval_sum = 0.0;

	memset( &m_estimated_rate, 0x0, sizeof( online_rate_entry_t));
	ret = 0;

#ifdef LWPR_TP
	//lwpr variable intialization
	memset(lwpr_x,0x0,sizeof(double)*23);
//	lwpr_x = (double*) malloc((nooffeatures)*sizeof(double)*23);
//	lwpr_y = (double*) malloc((double)*sizeof(double));
	lwpr_y = 0.0;
	//lwpr_y[0] = 0.0;
	lwpr_mse_T = 0.0;
	lwpr_mse_P = 0.0;
	lwpr_counter = 0;
	//c++ lwpr declarations..
	//lwpr_x(23);
	//lwpr_y(1);
//	lwpr_model(23,1);
//	lwpr_model.setInitD(50);
//	lwpr_model.setInitAlpha(250);
//	lwpr_model.wGen(0.2);
#endif
}

//MappleOnlineApplication destructor do nothing simply call the ns-log-function to print the message destructor
MappleOnlineApplication::~MappleOnlineApplication()
{
	NS_LOG_FUNCTION(this);
	fclose(fp);

#ifdef NORMALIZEDSAMPLES
	fclose(fp_error);
#endif
#ifdef LWPR_TP
	fclose(lwpr_file);
#endif
	free(oldSamples);
	free(recentSamples);
	free(testfeatures);
	//	fclose(fp1);
}
//might be I dont need this function, but at the moment I have it..
uint32_t MappleOnlineApplication::GetTotalRx(void) const
{
	NS_LOG_FUNCTION(this);
	return m_totalbytesRecv;
}
//on specific port who many neighbor nodes are conncected, and each of them have only this port to connect..
//even for this port we can count how many neighbors are connected..
Ptr<Socket> MappleOnlineApplication::GetListeningSocket (void) const
{
	NS_LOG_FUNCTION (this);
	return m_listningSocket;
}
//for how many nodes we accept their request for communication..
std::list<Ptr<Socket> > MappleOnlineApplication::GetAcceptedSockets (void) const
{
	NS_LOG_FUNCTION (this);
	return m_acceptedsocketList;
}
//simple IP4 address to ids mapping..
void MappleOnlineApplication::BuildIdMapp(void){

	m_totalNode = m_nodeContainer.GetN(); // get the number of nodes inside the node container
	//	NS_LOG_UNCOND("Number of nodes are:" << m_totalNode);
	Address source;
	uint32_t index;
	uint32_t ids;
	pair<Address, uint32_t > pair;
	for(index=0;index<m_totalNode; index++){
		m_actualNodes = m_nodeContainer.Get(index);
		source = m_actualNodes->GetDevice(0)->GetAddress();
		ids = m_actualNodes->GetId();
		pair.first = source;
		pair.second = ids;
		m_Mapp.insert(pair);
		//	std::cout<<"show at: "<<source <<" " << m_Mapp.at(source) << std::endl;
	}
}
//for DoDispose we not yet needed to destroy the allocated memory stuff.
//this is the main starting point of the application, even this StartApplication method is inherited from the
//object class, however Application class itself not provide the start() method to run the simulation..
//Secondly here we are going to deal with sockets because we need to send the summary packets and other related stuff..
//Even this would be valueable in future to send the learn model to other nodes using some Routing Protocol
void MappleOnlineApplication::StartApplication ()
{
	NS_LOG_FUNCTION (this);

	// at the moment we don't need this socket, we are just receiving the summaries through promiscuous mode..
	//At the first hand we need to creat the listning socket where we hear for incomming connections...
	if(!m_listningSocket){
		//basically GetNode() function will get the node to whom we will Start this application..
		//while m_tid would be the protocol type, like UDP, TCP etc
		m_listningSocket = Socket::CreateSocket(GetNode(), m_tid);
		//After creating socket we need to bind that socket on specific port with ordinary IPv4 address
		m_listningSocket->Bind (InetSocketAddress(Ipv4Address::GetAny(), 8000));
		//	m_listningSocket->Listen(); //now Listen() on created port All the time..this is blocking call
		//	m_listningSocket->ShutdownSend(); //We don't want to send something on this particular socket..only for listning
	}
	//at the moment am not interested to send anything..
	//however the socket assigned after the listening would be different see orignally BSD sockets..
	//because Listning itself blocking call so we have to assign some new socket for further communications ..
	//I really dont understand why we have SetRecvCallback() function before we Accepted the Socket..and Same function also called
	//from the HandleAccept() ..
	m_listningSocket->SetRecvCallback (MakeCallback (&MappleOnlineApplication::readSharedSamples, this));
	//this function have dual functionality and we need the second.. one where after listening the newly created socked created by
	//the second callback of the SetAcceptCallback() function...
	//		m_listningSocket->SetAcceptCallback (
	//					MakeNullCallback<bool, Ptr<Socket>, const Address &> (),
	//					MakeCallback (&MappleOnlineApplication::HandleAccept, this));
	//
	//		m_listningSocket->SetCloseCallbacks (
	//					MakeCallback (&MappleOnlineApplication::HandlePeerClose, this),
	//					MakeCallback (&MappleOnlineApplication::HandlePeerError, this));


	//at the moment we only interested to broadcast the summaries.
	//sendingSocket is used to send summaries packets
	if (!sendingSocket) {
		sendingSocket = ns3::Socket::CreateSocket(GetNode(),m_tid);
		//just using the different socket as compare to the application own socket..
		//InetSocketAddress remote = InetSocketAddress (Ipv4Address::GetBroadcast(), 81);
		destination = Ipv4Address::GetBroadcast();
		//				sendingSocket->Bind();
		//				sendingSocket->ShutdownRecv ();//We only interested to send data so shutdownrecv at the moment..
		sendingSocket->SetAllowBroadcast (true);
		//		sendingSocket->Connect(remote);
		//				sendingSocket->SetAttribute ("IpTtl", UintegerValue (1));	// Only on hop messages !!!

	}
	//creating socket for sending shared samples to other
	if(!sharedSocket){
		sharedSocket = ns3::Socket::CreateSocket(GetNode(),m_tid);
		//InetSocketAddress remote = InetSocketAddress (Ipv4Address::GetBroadcast(), 8000);
		//sharedSocket->Bind();
		sharedDestination = Ipv4Address::GetBroadcast();
		sharedSocket->SetAllowBroadcast(true);
		//sharedSocket->Connect(remote);
	}
	// Save a pointer to own mobility model
	//we need this.. whatever node who running this app, if we have movility model then we have to inform other nodes that
	//we are in moving state or standing state..etc.
	m_mobility =  GetNode ()->GetObject<RandomWaypointMobilityModel>();


	myid = GetNode()->GetId();
	snprintf(Nodeid, sizeof Nodeid, "%lu", (unsigned long)myid); /* Method 1 */
	saved_sample = "Sample_";
	saved_sample.append(Nodeid);
	saved_sample.append(".txt");
	printf("Saved_sample for node:%s\n",saved_sample.c_str());
	fp = fopen(saved_sample.c_str(),"w");

#ifdef NORMALIZEDSAMPLES //for training and prediction we need to give the model name
	model_name = "Model_";
	model_name.append(Nodeid);
	model_name.append(".model");
	printf("model name is:%s\n", model_name.c_str());
	read_model = "Model_";
	read_model.append(Nodeid);
	read_model.append(".model");
#endif
#ifdef TRAINING_PREDICTION
	fp_error = fopen("PredError.txt","w"); //writing the cumMSE for all the nodes
#endif

#ifdef LWPR_TP
	lwpr_init_model(&lwpr_model,23,1,"Mapple_lwpr");
	 /* Set initial distance metric to 50*(identity matrix) */
	lwpr_set_init_D_spherical(&lwpr_model,50);
	 /* Set init_alpha to 250 in all elements */
	   lwpr_set_init_alpha(&lwpr_model,250);

	   /* Set w_gen to 0.2 */
	   lwpr_model.w_gen = 0.2;

	lwpr_file = fopen("lwpr_error.txt", "w");
#endif
	//	predict_error = "PENode_";
	//	predict_error.append(Nodeid);
	//	predict_error.append(".txt");
	// Sava a pointer to MAC queue
	std::pair< Ptr< Ipv4 >, uint32_t > pair = interfaces.Get(GetNode()->GetId());
	Ptr< Ipv4 > ipv4 = pair.first;
	//smart pointer to wifinetdevice class
	Ptr< WifiNetDevice > wifiNetDevice = DynamicCast<WifiNetDevice>(ipv4->GetNetDevice(pair.second));	// get wifi net device
	Ptr< RegularWifiMac > regularWifiMAC = DynamicCast<RegularWifiMac>(wifiNetDevice->GetMac());		// get the MAC object
	macQueue = regularWifiMAC->GetDcaTxop()->GetQueue();// get the MAC queue object (tricky, required changing the RegularWifiMac class in order to make 'GetDcaTxop()' accessible...)
	//ok I get this point why we need mac queue.. even change GetDcaTxop() function private to public..IA
	// Also, send a promiscuous callback (when in unicast mode...)
	//if (SEND_UNICAST)//I dont know why Michal use this macro ..????
	BuildIdMapp();
	//ok we start receiving the traffic of local applications of this node that they want to send somewhere else..
	//wifiNetDevice->SetSendTrafficCallback(MakeCallback (&MappleOnlineApplication::HandleSendTraffic, this));
	wifiNetDevice->SetSendTrafficCallback(MakeCallback (&MappleOnlineApplication::HandleSendTraffic, this));

	wifiNetDevice->SetPromiscReceiveCallback(MakeCallback (&MappleOnlineApplication::HandlePromiscuous, this));
	//here might I need to use the schedulewithcontext() call
	//Simulator::Schedule (Seconds(2.0), &MappleOnlineApplication::SendingRate, this);

	//Simulator::ScheduleWithContext (GetNode ()->GetId (), Seconds (1.0), &MappleOnlineApplication::SendingRate);
	//	interval =  uniformvariable.GetValue(0.0,0.3);
	event_estimator = Simulator::Schedule(Seconds (1.0 + interval), &MappleOnlineApplication::RateEstimationTimer, this);
	summaryTimer = Simulator::Schedule (Seconds (1.0 + interval), &MappleOnlineApplication::sendSummaryToNeighbors,this);

	m_statusTimer = Simulator::Schedule (Seconds(NEIGHBOR_CHECK_INTERVAL+interval), &MappleOnlineApplication::NeighborStatusTimer, this);
	//Simulator::ScheduleNow(&MappleOnlineApplication::writeSamplesFile,this); //writing into sample file);
}
//nothing reading at the moment// even we might not need this function..
void MappleOnlineApplication::readSharedSamples (Ptr<Socket> socket)
{
	NS_LOG_FUNCTION (this << socket);
#ifdef SHARED_SAMPLES
	Ptr<Packet> packet;
	Address fromAddress;
	//	while (packet = socket->RecvFrom (fromAddress))
	//	{
	packet = socket->RecvFrom (fromAddress);
	//	printf("The packet size is:%d,sizof:%d\n",packet->GetSize(),sizeof(packet));
	if (packet->GetSize () == 0)
	{ //EOF
		printf("Pakcet for shared sample is 0\n");
		//		break;
	}

	if (InetSocketAddress::IsMatchingType (fromAddress))
	{
		//m_totalRx += packet->GetSize ();
		Ipv4Address address = InetSocketAddress::ConvertFrom (fromAddress).GetIpv4();
		NS_LOG_LOGIC (GetNode()->GetId()  << " (" << interfaces.GetAddress(GetNode()->GetId()) << "): Received " << packet->GetSize () << " bytes from " << " [" << address << "]");


		uint8_t buf[packet->GetSize()];
		packet->CopyData(buf,packet->GetSize());
		shared_samples_t *recv_sharedsamples = (shared_samples_t*) buf;
		//printf("Sizeof() packet:%d------------------------\n",sizeof(packet));
		//memset(&packet,0x0,sizeof(packet));
		if(recv_sharedsamples->msg_type == SHAREDSAMPLES){
			if(recv_sharedsamples->node_id != myid){
				ss_counter += 1;
				memcpy(ssamples.samples_array,recv_sharedsamples->sharedfeatures,sizeof(double)*24);

				samples_queue.push(ssamples);
				//						samples_queue.push(recv_sharedsamples->sharedfeatures);
				//						printf("Ret: %d recve Shared samples from :%u count:%u\n",ret,recv_sharedsamples->node_id,ss_counter);

//				printf("Sample Received:\n");
//				for(int i=0;i<=22;i++){
//					printf("%f ",recv_sharedsamples->sharedfeatures[i]);
//				}
//				printf("\n");
				//					}
				//					memcpy(comparing,recv_sharedsamples->sharedfeatures,sizeof(double)*23);
				//					}
				memset(recv_sharedsamples,0x0,sizeof(shared_samples_t));
			}
		}
	}

	//	}
#endif

}

void MappleOnlineApplication::HandleAccept (Ptr<Socket> s, const Address& from)
{
	NS_LOG_FUNCTION (this << s << from);
	s->SetRecvCallback (MakeCallback (&MappleOnlineApplication::readSharedSamples, this));
	m_acceptedsocketList.push_back(s);
}

void MappleOnlineApplication::HandlePeerClose (Ptr<Socket> socket)
{
	NS_LOG_INFO ("MappleOnlineApplication, peerClose");
}

void MappleOnlineApplication::HandlePeerError (Ptr<Socket> socket)
{
	NS_LOG_INFO ("MappleOnlineApplication, peerError");
}
//Note// we might not need additional arguments for this function because we only concern with packet size, and we might be
//finally through the size of the packet from wifi-net-device.cc instead we reference to the packet from the lower layer, which
//could be dangeroud becuase what if packet over write or something else could happen..
bool MappleOnlineApplication::HandleSendTraffic(const Address & dest, uint32_t lenth, double sent_time, uint16_t protocol) {
	//NS_LOG_UNCOND("RECEIVING TRAFFIC INSIDE THE HANDLESENDTRAFFIC FUNC\n");
	//float currentRate;
	m_estimated_rate.bytes_sent += lenth;
	sent_time_counter += sent_time;
	NS_LOG_INFO("Sent_Time is: " << sent_time);
	if(hasEMA == 0){
		last_sent_time = sent_time;
		hasEMA = 1;
		m_estimated_rate.packet_len = lenth;
		//		flag = false;
	}else {
		float seconds = sent_time - last_sent_time; // / sim_ticks_per_sec();
		//	NS_LOG_LOGIC("\nSeconds calculated: " << seconds);

		em_interval_sum += seconds;
		// EMA_RATE_W = 0.013
		m_estimated_rate.packet_len = (1.0 - EMA_RATE_W) * m_estimated_rate.packet_len + (EMA_RATE_W)*lenth;
		NS_LOG_LOGIC("\nRateEstimator Estimate the packet len: " << m_estimated_rate.packet_len);


		if( hasEMA == 1)
		{
			m_estimated_rate.interval = seconds;
			m_estimated_rate.rate =(int) ceil(m_estimated_rate.packet_len / m_estimated_rate.interval);

			hasEMA = 2;
		}else	//after we receive third packet from wifi-netdevice.cc this will always execute
		{
			m_estimated_rate.interval = ((1.0 - EMA_RATE_W) * m_estimated_rate.interval) + ((1.0 * EMA_RATE_W) * seconds);

			m_estimated_rate.rate = (int) ceil(m_estimated_rate.packet_len / m_estimated_rate.interval);
			NS_LOG_LOGIC("\nEstimated rate is: " << m_estimated_rate.rate);

		}
		last_sent_time = sent_time;
		///need to work on this logic///
		//	flag = false;
		//if(sent_time_counter == Simulator::Now().GetSeconds()){
		//	Simulator::ScheduleNow(&MappleOnlineApplication::RateEstimationTimer);
		//	}
	}
	sent_traffic_counter += 1;
	//printf("sent traffic counter:%u node:%u at time:%f \n",sent_traffic_counter, myid,Simulator::Now().GetSeconds());
	return true;
}
//this is promisscuous call just change the name from the michal function
bool MappleOnlineApplication::HandlePromiscuous(Ptr<NetDevice> device, Ptr<const Packet> packet, uint16_t protocol,
		const Address & sender, const Address & receiver, enum NetDevice::PacketType packetType) {
	NS_LOG_FUNCTION_NOARGS ();
	//	printf("PACKET RECEIVED IN PROMISSCUOUS FUNCTION\n");


	// Here, we deal only with packets addressed to other hosts
	//	if (packetType != NetDevice::PACKET_OTHERHOST) //Otherhost mean this packet not belong to us..
	//	return true;				//well for online we need this packet either its belong to us or not.. see mapple paper
	// We are only interested in Ipv4 packets, and not the packets from lower layers...
	if (protocol != Ipv4L3Protocol::PROT_NUMBER) //we need this check because the lower like mac have many packets their own, but we need only which have some payload or packet belong to application layer
		return true;

	//	NS_LOG_UNCOND("PACKET RECEIVED IN PROMISSCUOUS FUNCTION Again");

	Ipv4Header header;
	packet->PeekHeader(header); //not remove the header.

	//Ipv4Address dst = header.GetDestination ();
	//Ipv4Address origin = header.GetSource ();
	uint32_t id=0;
	id = m_Mapp.at(sender); //here we have complete map of the Mac48Address to ids..
	//here I might to perform the check inorder to escape from those summaries which I am currently broadcasting for others
	//so simpley ignor those packets which are sending by this node.. Ask Michal?
	Ptr<Packet> p2 = packet->Copy();	// need a copy to remove headers (can not operate on 'packet' as it is declared constant...)
	Ipv4Header ipHeader;
	p2->RemoveHeader(ipHeader);			// remove IP header
	LlcSnapHeader macHeader;
	p2->RemoveHeader(macHeader);		// remove MAC layer header

	uint32_t s = p2->GetSize();			// get data, finally
	//printf("packetsize is:%u \n",s);
	uint8_t buf[s];
	p2->CopyData(buf,p2->GetSize());
	summary_msg_t* summary = (summary_msg_t*) buf;
	//
	//	if(p2->GetSize() > 130 ){
	//
	//	}else{
	//		uint32_t *countp;
	//		uint32_t s2 = p2->GetSize();
	//		uint8_t buf2[s2];
	//		p2->CopyData(buf2,p2->GetSize());
	//		countp = ((uint32_t*)(buf2));
	//
	//		printf(" ----------------NODE : %u Sendpacket sequence: %u  At time :%f and Received At: %u-----------------\n ", id ,*countp , Simulator::Now().GetSeconds() ,  GetNode()->GetId());
	//
	//	}


	// Get RSSI of received message
	double rssi = 0;
	double snr = 0;
	PhySignalTag pst;
	if(packet->PeekPacketTag(pst)){		// the RSSI value is stored in a packet tag added my (modified) 'YansWifiPhy::EndReceive' function
		rssi = pst.GetPwr();
		snr = pst.GetSnr();	// As a 'bonus', also the SNR value is obtained here...
	} else {
		cout << "A serious error - no packet tag with the RSSI value available !!!"<< endl;
	}
	////////////////////////////////As Michal either we need this or not/////////////////////////////////////
	//		double rat = 0;
	//		double ht = 0;
	//	MacDelaysTag mdt;
	//	if(packet->PeekPacketTag(mdt)){		// the RSSI value is stored in a packet tag added my (modified) 'YansWifiPhy::EndReceive' function
	//		rat = mdt.GetRequestAccessTime();
	//		ht = mdt.GetHandshakeTime();	// As a 'bonus', also the SNR value is obtained here...
	//	} else {
	//		cout << "A serious error - no packet tag with MAC delays values available !!!"<< endl;
	//	}
	//		uint32_t qs = 0;
	//	MacQueueTag mqt;
	//	if(packet->PeekPacketTag(mqt)){		//
	//		qs = mqt.GetQueueSize();
	//	} else {
	//		cout << "A serious error - no packet tag with MAC queue values available !!!"<< endl;
	//	}
	handlePacketReception(summary,id,rssi,snr);
	////////////////////////////////As Michal either we need above block or not/////////////////////////////////////

	return true;
}
//1. call from the handlePromisscuous()
void MappleOnlineApplication::handlePacketReception(summary_msg_t* recv_summary, uint32_t id, double rssi, double snr) {

	//uint32_t myid = GetNode()->GetId();

	//	double recv_time = Simulator::Now().GetSeconds();
	bool envChanged = false;

	Neighbors[id].last_message_time = Simulator::Now().GetSeconds();;
	Neighbors[id].last_rssi = rssi;
	Neighbors[id].last_snr = snr;

	//see at the end of program we again setting the same value..
	if(!isNeighborStatusActive(id)){ //
		envChanged = true; //if it is not active neighbor then environment change at the receiver, because we modified the nieghbor table
		//		printf("\n Neigbor:%u status is not active at:%u\n",id,myid);
		setNeighborStatus(id, STATE_ACTIVE); //for this particular neighbor we will set state_active
	}else if(recv_summary->msg_type == SUMMARY_PACKET){
		//		printf("We Received the summary packet\n");
		if( myNeighborsRateChanged(recv_summary, id)){
			envChanged = true;
			//			printf("\nmyNeighbor:%u Ratechanged at:%u\n",id,myid);
		}

	}else{
		NS_LOG_LOGIC("Received Normal Packet"<<id<<"\n");
	}
	// If my environment changed - close all samples
	//page 16: local environment change when
	//1. Estimated sending rate change
	//2. its Neighbor table Modified ( enter some new node in it)
	//3. Sending rate of the neighbors changed
	// for the local movements and neighbor movements we are not considering the environment change yet!!!!!!
	if( envChanged)
	{
		//	printf(" ENV CHANGE from handle packet reception.....\n");
		environmentChange(Simulator::Now().GetSeconds());
	}

	if( recv_summary->msg_type == SUMMARY_PACKET)
	{
		//		printf("\n summary Recv at:%u ontime:%lf  FROM: %u \n",myid,Simulator::Now().GetSeconds(),id);
		//reportSummary(recv_summary, from, recv_time); //show the messages what we received
		processNeighborSummary(recv_summary,id,rssi,snr); //it seems to be awkward to passing so many time the same parameter, but sake of clearity
		//memset(recv_summary,0,sizeof(summary_msg_t));
	}else
	{
		// count it : if its other than they the summary packet, and according to flow chart of the protocol, we are not allow to
		//start sampling the neighbor node. At the reception of the summary packet we start sampling...However this assumption is
		// great that, all the srrounding nighbors running the mapple protocol. but what if neighbor nodes dont have the mapple, can we still consider them. reference for future?
		if( hasOpenSample(id))
		{
			Samples[id].received += 1;
			Samples[id].sum_rssi += rssi;
			//printf("Sum_rssi is:%d Rcounter:%u and rssi:%f\n",Samples[id].sum_rssi,Samples[id].received,rssi);
			Samples[id].sum_snr += snr;
			//	Samples[id].last_rssi = rssi; //this is for calculating the rssi rates of the particular neighbors
		}
	}
	//Eduardo Idea to place them here, but I am getting them on top!
	//	//even if above condition become false, we still count for it, and next time when we will open sample with the help of summary packet, these values would count
	//	Neighbors[id].last_message_time = recv_time;
	//	Neighbors[id].last_rssi = rssi;
	//	Neighbors[id].last_snr = snr;


}

//1. call from handlePacketReception()
void MappleOnlineApplication::processNeighborSummary(summary_msg_t* recv_summary, uint32_t id, double rssi, double snr)
{
	//uint32_t myid = GetNode()->GetId();
	//printf("Processing the Summary\n");
	// First sample for neighbor "id"
	if(isFirstPacketSample(id))
	{
		openSamplesDatabase(recv_summary,id,rssi,snr);

		//		printf("\n Open Sample Done at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
	}else
	{
		// Summaries can tell if a current sample must be closed
		if( isNeighborEnvironmentChanged(recv_summary, id)){ //Actually we are performing 3 check in this if. status, neighbors sending rates, and last_update time
			//			printf("\n if isNebENC at:%u at time:%lf for node:%u \n" ,myid,Simulator::Now().GetSeconds(),id);
			//if above check true then close all the samples, and open sample for this message
			closeSample(id);
			//			printf("\n Opening after closing at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
			openSamplesDatabase(recv_summary,id,rssi,snr); //this function call at first time when we received the first sample after setting the active sate of the neighbors
		}else{ //we simple increment the sample counts
			//			printf("\n\n updatesamples at:%u at time:%lf for node:%u \n" ,myid,Simulator::Now().GetSeconds(),id);
			updateSamplesDatabase(recv_summary,id,rssi,snr);// this for second time to receive the samples and onward..

		}
	}
	//printf("proccessedSummary received current rate:%d " ,recv_summary->my_current_rate);
	Neighbors[id].neighbors_rate = recv_summary->my_current_rate;
	Neighbors[id].last_update = recv_summary->last_update; // help to discard last sample, see orignal idea of Mapple, to drop last message.
	//because if this one is the new sample then, we have check that if last update time is less then the received, then its means that sender
	//have a change in its environments..etc
}

//1. call from handlePacketReception() 2 time
//2. call from closeAllSamples()
bool MappleOnlineApplication::hasOpenSample(uint32_t id)
{
	return (Samples[id].open == true);
}
//1. call from processNeighborSummary() //how ever Michal doing implementation for openSamplesDatabase and updatSample within one function

void MappleOnlineApplication::openSamplesDatabase(summary_msg_t* recv_summary, uint32_t id, double rssi, double snr)
{
	uint32_t i;
	//double current_time = Simulator::Now().GetSeconds();

	// Consider it as a normal packet
	Samples[id].last_rssi = Neighbors[id].last_rssi;
	Samples[id].received = 1;
	Samples[id].sum_rssi = rssi;
	Samples[id].closed_sum_rssi = 0;
	Samples[id].sum_snr = 0;
	Samples[id].closed_sum_snr = 0;
	Samples[id].open = true;
	for(i=0; i<5; i++){
		Samples[id].sender_status[i] = recv_summary->source_status[i]; //how many active neighbors sender hold
		Samples[id].receiver_status[i] = my_status[i]; //how many neighbors we hold
		//this was the orignal idea of eduardo to update receiver status only at sample opening, and not update during updatesamples
		//in contrary Michal updating both places.. I think this will not harm anything
	}
	Samples[id].activeNeighbors = recv_summary->activeNeighbors; //number of active nighbors of the source//seems redundant ...
	//here we only interested to take care number of neighbors of particular node, but not which one.
	for(i=0; i<recv_summary->activeNeighbors; i++){
		Samples[id].sender_rssi_data[i] = recv_summary->neighbors_rssi_data[i];
		Samples[id].sender_rate_data[i] = recv_summary->neighbors_rate_data[i];
		Samples[id].sender_snr_data[i] = recv_summary->neighbors_snr_data[i];
	}
	//michal additional features.. They are seems to be good features for defining the routing paths..
	//	Samples[id].avgTravelTime = (current_time - recv_summary->sendTime);
	//	Samples[id].queueSize = qs;
	//	Samples[id].avgSendTime = (current_time - rat);	// calculated starting from the moment when dca-txop requested for channel access for this packet...
	//	Samples[id].avgSendTimeWithoutBackoff = (current_time - ht);	// calculated starting from the moment when packet was dequeued from the MAC queue


	Samples[id].sender_rate = recv_summary->my_current_rate;
	Neighbors[id].last_update = recv_summary->last_update; //after environment change we set this last_update.
	Samples[id].first_counter = recv_summary->sent_traffic_counter; //for first counter and sent_traffic_counter_lastclose at opening of sample we have sent_traffic_counter for both
	Samples[id].sent_traffic_counter_lastclose = recv_summary->sent_traffic_counter;
	Samples[id].open_time = Simulator::Now().GetSeconds();

	NS_LOG_INFO("LinkQualityEstimator "<< Samples[id].open_time << " Opening sample for " << id << " with last_update " << recv_summary->last_update);

}

//1. call form processNeighborSummary() // we can deal first and later sample in one function as Michal doing.!
void MappleOnlineApplication::updateSamplesDatabase(summary_msg_t* recv_summary, uint32_t id, double rssi, double snr)
{
	uint32_t i, packets_sent;
	//double current_time = Simulator::Now().GetSeconds();
	//	printf("-++++++++++++++++++++++++++Update Sample+++++++++++++++++++++++++++++++++++++++++\n");
	Samples[id].last_rssi = Neighbors[id].last_rssi; //help to calculate the neighbors rssi rate calculations..
	// Get updated rssi data
	for(i=0; i< recv_summary->activeNeighbors; i++){
		Samples[id].sender_rssi_data[i] = recv_summary->neighbors_rssi_data[i];
		Samples[id].sender_snr_data[i] = recv_summary->neighbors_snr_data[i];
		//I do beleive that snr is important feature, mixture of snr and rssi could decide the best link quality..
		Samples[id].sender_rate_data[i] = recv_summary->neighbors_rate_data[i];
	}
	//since sent_traffic_counter is == to number of packets applications sent out at the sender side, which is available us through
	//summary packet.. and the sample.sent_traffic_counter_lastclose gives us the pervious sent_traffic_counter.. Their minus provide the total packet sent
	//by the sender node

	packets_sent = recv_summary->sent_traffic_counter - Samples[id].sent_traffic_counter_lastclose;
	//printf("UDS->sent_traffic_counter:%u - lastclos:%u from:%u at me:%u time:%f and packet_sent:%u\n",recv_summary->sent_traffic_counter,Samples[id].sent_traffic_counter_lastclose,id,myid,Simulator::Now().GetSeconds(),packets_sent);
	NS_LOG_INFO("LinkQualityEstimator " << Simulator::Now().GetSeconds()<<" Updating sample. Closing " <<Samples[id].received<<" packets received and " <<	packets_sent<< " packets sent\n");

	//here we again updating the last close with sender's sent_traffic_counter
	Samples[id].sent_traffic_counter_lastclose = recv_summary->sent_traffic_counter;
	Samples[id].last_close_time = Simulator::Now().GetSeconds();
	Samples[id].packets_recv_closed += Samples[id].received; //in handlePacketReception() if samples are open for id then we increament the received by 1
	Samples[id].packets_sent_closed += packets_sent;

	//	printf("UDS: packet_recv_closed:%u, sent_closed:%u for:%u atme:%u time:%f\n",Samples[id].packets_recv_closed,Samples[id].packets_sent_closed,id,myid,Simulator::Now().GetSeconds());


	Samples[id].closed_sum_rssi += Samples[id].sum_rssi; //when receiving the last packe while closing the sample we have the opertunity to store them..etc
	//	printf("Sample.ClosedSumrssi is:%d and SampleReceived:%u\n",Samples[id].closed_sum_rssi,Samples[id].packets_recv_closed);
	Samples[id].closed_sum_snr += Samples[id].sum_snr;

	//at the moment we dont need this fucntionality..
	//	Samples[id].avgTravelTime += (current_time - recv_summary->sendTime);
	//	Samples[id].queueSize += qs;
	//	Samples[id].avgSendTime += (current_time - rat);	// calculated starting from the moment when dca-txop requested for channel access for this packet...
	//	Samples[id].avgSendTimeWithoutBackoff += (current_time - ht);	// calculated starting from the moment when packet was dequeued from the MAC queue

	for(i=0; i<5; i++){
		Samples[id].sender_status[i] = recv_summary->source_status[i]; //how many active neighbors sender hold
		Samples[id].receiver_status[i] = my_status[i]; //how many neighbors we hold at particular sender neighbor id
	}
	//	memcpy(Samples[id].sender_status,recv_summary->source_status,sizeof(uint32_t)*5);
	//	//	memcpy(Samples[id].receiver_status, &my_status[0],sizeof(uint32_t)*5); //Michal did this, but seems to be awkward, that
	//	//we reassigned the my_status value through indexes.
	//
	//		Samples[id].receiver_status[0] = my_status[0];
	//		Samples[id].receiver_status[1] = my_status[1];
	//		Samples[id].receiver_status[2] = my_status[2];
	//		Samples[id].receiver_status[3] = my_status[3];
	//		Samples[id].receiver_status[4] = my_status[4];

	if(Samples[id].packets_sent_closed > MAX_PACKETS_CLOSED)
	{
		closeSample(id);
		NS_LOG_LOGIC("LinkQualityEstimator " << Simulator::Now().GetSeconds() << " Sample has MAX_PACKETS_CLOSED\n");
	}

	// Start again counting
	Samples[id].received = 1;
	Samples[id].sum_rssi = rssi;
	Samples[id].sum_snr = snr;
}



//1. call from processNeighborSummary()
bool MappleOnlineApplication::isNeighborEnvironmentChanged(summary_msg_t *recv_summary, uint32_t id)
{
	uint32_t i;
	//uint32_t myid=GetNode()->GetId();

	NS_LOG_LOGIC("LinkQualityEstimator"<< Simulator::Now().GetSeconds() << " Checking neighbor environment with summary\n");
	for(i=0; i< 5; i++)
	{
		if( Samples[id].sender_status[i] != recv_summary->source_status[i]){
			//			printf("\nSender status not match at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
			return true;
		}
	}
	for(i=0; i< recv_summary->activeNeighbors; i++)
	{

		if( Samples[id].sender_rate_data[i] != recv_summary->neighbors_rate_data[i]){
			//			printf("\nSenders rates not match at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
			return true;
		}
	}
	if(getCurrentRssiRate(id) !=  getOldRssiRate(id)){
		//		printf("\nSenders RSSI not match at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
		return true;
	}
	// Neighbor environment changed since last
	// received summary.last_update time is currently set so its greater than the last_update time
	//2.99 < 2.99 condition false and exit this check
	if( Neighbors[id].last_update < recv_summary->last_update){
		//		printf("\nNeb last update not match at:%u at time:%lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),id);
		return true;
	}
	return false;
}
//environmentChange() mean close all the sample..CHECK the page 16 description of the MAPPLE protocol...
//1. call from rateEstimationTimer()
//2. call from handlePacketReception()
//3. call from NeighborStatusTimer()
void MappleOnlineApplication::environmentChange(double t)
{
	NS_LOG_INFO("ENVIRONMENT CHANGE at time:" << t <<"- Closing all samples at node:" << GetNode()->GetId() <<"\n");
	//NS_LOG_UNCOND("Environment Changed at "<< t <<" Forcing SUMMARY send\n");
	//Simulator::Cancel(summaryTimer);
	closeAllSamples();
	last_update = t;
	if (IMMEDIATE_SUMMARY){
		//call SummaryTimer.startOneShot(1 );
	}
	// as Eduardo did for summay timer to fire on request.. I also need the same time in ns3..
}
//1. call from processNeighborSummary();
//2. call from closeAllSample()
void MappleOnlineApplication::closeSample(uint32_t id)
{
	NS_LOG_LOGIC("close sampe at time: " << Simulator::Now().GetSeconds() << " for %u\n" << id);
	
	if( Samples[id].packets_sent_closed >= MIN_PACKETS_CLOSED)
	{
		//here we have trade off, hum much we want to gain, if we received lower than the 50 samples we simple ignore all the previous
		//sample of this id..

		printSample(id);
		NS_LOG_INFO("PRINTING SAMPLES");
	}

	//inorder to complete the partial samples which were less than the MIN_PACKET_CLOSED, so this is the best place to start keep tracking them
	//Instead we wash out particular sample of the sender we will store them and if in future we encounter the similar environment then we will start
	//over where we left and once we will cross the MIN_PACKET_CLOSED threshold we will use that sample for learning.
	//Furthermore we need seprate data structure to handle them, and we also want their discrete values plus some aditional parameter like how many
	//of them we received so far.

	memset(&Samples[id], 0x0, sizeof(online_sample_t));
	/* Make sure open flag is set to false */
	Samples[id].open = false;
	Samples[id].packets_sent_closed = 0;
	Samples[id].packets_recv_closed = 0;
	Samples[id].closed_sum_rssi = 0;
	Samples[id].sum_rssi = 0;
	Samples[id].sum_snr = 0;
	Samples[id].closed_sum_snr = 0;
	Samples[id].last_rssi = 0;
	Samples[id].open_time = 0;
	Samples[id].last_close_time = 0;
	Samples[id].sent_traffic_counter_lastclose = 0;
}

//1. call from environmentChange()
void MappleOnlineApplication::closeAllSamples(void)
{
	uint32_t i =0;
	for(i=0;i < MAX_NODE_ID; i++)
	{
		if( hasOpenSample(i))
		{
			closeSample(i);
		}
	}
}
//this function since checking the nodes own neighbor rates, if they have change we should close the samples
//1. call from HandlePromiscuous()
bool MappleOnlineApplication::myNeighborsRateChanged(summary_msg_t *recv_summary, uint32_t id)
{

	if( recv_summary->my_current_rate != Neighbors[id].neighbors_rate)
	{
		//		printf("\nSender rate is different then the previous one\n");
		return true;
	}
	return false;
}
//1. call from HandlePacketReception()
//2. call from NeighborStatusTimer.fired
bool MappleOnlineApplication::setNeighborStatus(uint32_t id, uint8_t state)
{
	if( state == STATE_NOT_ACTIVE)
	{
		Neighbors[id].state = STATE_NOT_ACTIVE;
		clear_bit(my_status[id/32],id % 32);
	}else if(state == STATE_ACTIVE)
	{
		Neighbors[id].state = STATE_ACTIVE;
		set_bit(my_status[id/32],id % 32); //what that mean?
	}
	return true;
}
//check that, for the sender, is this particular sender (neighbor) active or not
bool MappleOnlineApplication::isNeighborStatusActive(uint32_t id)
{

	return (Neighbors[id].state == STATE_ACTIVE);
}
//check that, is it first packet from the sender(neighbor) or not
bool MappleOnlineApplication::isFirstPacketSample(uint32_t id)
{
	return (Samples[id].open == false);
}
//Eduardo Implemention for RateEstimationTimer, its independent function, and no one call it!
void MappleOnlineApplication::RateEstimationTimer(void){

	float factor = 1024.0 / RATE_ESTIMATION_INTERVAL ; //RATE_ESTIMATION_INTERVAL = 4096, // milliseconds ?? 0.250
	float current_rate;
	oldRate = getCurrentSendingRate();
	//	printf("OldRate -< Getcurrentrate: %d \n", getCurrentSendingRate());
	oldBytesSent = m_estimated_rate.bytes_sent;
	oldEstRate = m_estimated_rate.rate;

	current_rate = factor * m_estimated_rate.bytes_sent; //125 for rate 1, we sent 500 bytes/second

	if (hasEWMA == 0){
		m_estimated_rate.rate = (int) ceil(current_rate);
		//printf("Ratetimerfired First time and rate is: %u at nod3:%u\n",m_estimated_rate.rate,myid);
		hasEWMA = 1;
		//NS_LOG_INFO("\nHasEWMA: " << hasEWMA);
		//	Simulator::Schedule (Seconds(1.0), &MappleOnlineApplication::RateEstimationTimer);

	}else{
		m_estimated_rate.rate = (int) ceil( ((1.0 - EWMA_RATE_W) * m_estimated_rate.rate) + (1.0*(EWMA_RATE_W)*current_rate));
		//	printf("Ratetimerfired and rate is: %u at node:%u \n",m_estimated_rate.rate,myid);
		//NS_LOG_INFO("\nHasEWMA: " << hasEWMA);
	}
	m_estimated_rate.bytes_sent = 0;
	NS_LOG_LOGIC("\nRate Estimator at time: " << Simulator::Now().GetSeconds() << " OldBytesSent: " << oldBytesSent << " OldEstimated rate: "<<oldEstRate
			<< " New ByteRate: \n" << m_estimated_rate.rate	);

	if( oldRate != getCurrentSendingRate())
	{
		printf("At time %lf OldRate %d  is Different than current Sending rate %d at Node:%u : \n",  Simulator::Now().GetSeconds(), oldRate,  getCurrentSendingRate(),GetNode()->GetId());
		environmentChange(Simulator::Now().GetSeconds());
	}
	if(Simulator::Now().GetSeconds() < m_Simulationtime){
		//for rate estimation we do not need randomness
		//		interval =  uniformvariable.GetValue();
		event_estimator=	Simulator::Schedule (Seconds(1.0+interval), &MappleOnlineApplication::RateEstimationTimer,this);
		//Sending the summary packet every time when we estimating our application sending rates, but here I need oneshot
		//timer as Eduardo doing in his implementation..becasue we need to fire the summaries while our neighbors rate chages etc..
		//		Simulator::Schedule (Seconds(1.0), &MappleOnlineApplication::sendSummaryToNeighbors, sendingSocket);
	}
}
//findind the RSSI rates, we need two rates old and current, these two function seems to
//redundant to me.. as compare to getSendersNeighborsRssiRates(rssi), here we just pass
// rssi value find its rate, I will decide later how to deal with them?
//1. call from isNeighborEnvironmentChanged()
uint8_t MappleOnlineApplication::getOldRssiRate(uint32_t id){

	if( Samples[id].last_rssi > RSSI_RATE_LIMIT_MED){
		return RSSI_RATE_HIGH;
	}else if(Samples[id].last_rssi <= RSSI_RATE_LIMIT_MED && Samples[id].last_rssi > RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_MED;
	}else if(Samples[id].last_rssi < RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_LOW;
	}

	return 0;
}
//1. call from isNeighborEnvironmentChanged()
uint8_t MappleOnlineApplication::getCurrentRssiRate(uint32_t id){

	if( Neighbors[id].last_rssi > RSSI_RATE_LIMIT_MED){
		return RSSI_RATE_HIGH;
	}else if(Neighbors[id].last_rssi <= RSSI_RATE_LIMIT_MED && Neighbors[id].last_rssi > RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_MED;
	}else if(Neighbors[id].last_rssi < RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_LOW;
	}

	return 0;
}
//1. call from RateEstimationTimer()
//2. call from sendSummaryToNeighbors()
//3. call from printSamples()
uint8_t MappleOnlineApplication::getCurrentSendingRate(void){

	if( m_estimated_rate.rate < LQE_RATE_LIMIT_MED){
		return LQE_RATE_LOW;
	}else if( m_estimated_rate.rate >= LQE_RATE_LIMIT_MED && m_estimated_rate.rate < LQE_RATE_LIMIT_HIGH){
		return LQE_RATE_MED;
	}else if( m_estimated_rate.rate > LQE_RATE_LIMIT_HIGH){
		return LQE_RATE_HIGH;
	}

	return 0;
}
//filling the summary packet and send them//
//1. schedual from Start Application
//2. schedual itself as in loop
void MappleOnlineApplication::sendSummaryToNeighbors(void){

	NS_LOG_FUNCTION_NOARGS ();
	//printf("We are in sendSummaryToNeighbors()\n");
	int j,k;
	//NS_LOG_UNCOND ("Sending packet From EQSP at: " << Simulator::Now().GetSeconds());
	//	NS_ASSERT (m_sendEvent.IsExpired ());

	summary_msg_t send_summary;
	//send_summary = 0;
	//printf("Summary structure initialized\n");
	memset(&send_summary,0x0,sizeof(summary_msg_t));
	//	send_summary->sender_id = GetNode()->GetId(); //actuall we dont need this, we already build macaddress to id map..
	send_summary.msg_type = SUMMARY_PACKET;
	//	printf("Summary type setting done\n");
	send_summary.counter = ++counter; // this is just summary packet counter
	for(int j=0; j< 5;j++){
		send_summary.source_status[j] = my_status[j];
	}
	send_summary.last_update = last_update; //this time value come from environmentChange trigger..
	//everytime we have to send this last_update time, unless environment trigger and we update this value..
	//suppose we have 2.999 value

	send_summary.my_current_rate = getCurrentSendingRate(); //my own applications sending rates..
	//	printf("Summary filled with current rate:%d " ,send_summary.my_current_rate);

	k=0;
	for(j=0;j < MAX_NODE_ID; j++)
	{
		if(isNeighborStatusActive(j))
		{
			// after every seconds we are sending the average rssi of the particular neighbor id
			send_summary.neighbors_rssi_data[k] =getNeighborRSSI(j);// (Samples[j].sum_rssi) / (Samples[j].received)+1  ; //
			send_summary.neighbors_snr_data[k] = getNeighborSNR(j);//(Samples[j].sum_snr) / (Samples[j].received)+1 ;//
			send_summary.neighbors_rate_data[k] = Neighbors[j].neighbors_rate;
			//	printf("SUMMARY RSSI and RATE DATA for :%d : Rssi: %d, Snr:%d at Rate: %d \n", j , send_summary.neighbors_rssi_data[k], send_summary.neighbors_snr_data[k], send_summary.neighbors_rate_data[k]);
			k+=1;
		}
		if( k >= MAX_REPORTED_NEIGHBORS){
			NS_LOG_UNCOND("ERRORS"<< Simulator::Now().GetSeconds() <<" ERROR - Ups, too many neighbors\n");
		}
	}
	send_summary.activeNeighbors = k; // this parameter is interesing one! why because we are only interested in numbers of neighbors, but
	//not in particular? but could be criticle if we take care of individual one?
	//NS_LOG_UNCOND("SUMMARY activeNeighbors " << send_summary.activeNeighbors);
	//	printf("\nSumary Send:Active Neighbors At Node:%u: %d\n", GetNode()->GetId() ,k);
	send_summary.sent_traffic_counter = sent_traffic_counter; //sent_traffic_counter is all those numbers of packets which we sent from begining of time
	NS_LOG_INFO("SUMMARY counter "<< send_summary.counter <<" sent_traffic_counter :" << send_summary.sent_traffic_counter);

	//affter filling the packet we need to send the packet.. I will see this later!
	//printf("Summary is filld and ready to send\n");

	Ptr<Packet> packet = Create<Packet> (reinterpret_cast<const uint8_t*>(&send_summary),sizeof(summary_msg_t));
	//destination = Ipv4Address::GetBroadcast();
	sendingSocket -> SendTo(packet,0,InetSocketAddress(destination,portBase)); //
	//	printf("Enqueue summary send done");

	if(Simulator::Now().GetSeconds() < m_Simulationtime){

		summaryTimer = Simulator::Schedule (Seconds(1.0+interval), &MappleOnlineApplication::sendSummaryToNeighbors,this);
	}
}
//At the moment same implementation as the Eduardo did.! After every 2.048 seconds we will check either our neighbor list is same
//or changed..
void MappleOnlineApplication::NeighborStatusTimer(void){
	int i;
	bool envChange = false;
	//uint32_t myid = GetNode()->GetId();

	//NS_LOG_INFO("LinkQualityEstimator  Checking Neighbors ");
	for(i=0; i< MAX_NODE_ID; i++)
	{
		if( Neighbors[i].state == STATE_ACTIVE)
		{
			//max neighbor time to live is 2 seconds.. If we are not receiving anything from neighbor "i" we will
			//remove that neighbor after 2 seconds..see the orignal description from the paper!
			//however if my randomwaypoint modility speed increase, then this condition trigger allot even I increase ttl 3 but
			//average difference I found is 6 seconds I dont know why? ask Michal?
			if( (Simulator::Now().GetSeconds() - Neighbors[i].last_message_time + interval) > MAX_NEIGHBOR_TTL){
				//	NS_LOG_INFO("LinkQualityEstimator "<< Simulator::Now().GetSeconds() <<"  Removing inactive neighbor " << i);
				//	printf("\nNebStatusTimer fired at:%u at time:%lf - %lf for node:%u\n" ,myid,Simulator::Now().GetSeconds(),Neighbors[i].last_message_time,i);
				setNeighborStatus(i, STATE_NOT_ACTIVE);
				envChange = true;

			}
		}
	}
	//if neighbor state is change then we must trigger the environmentChange()
	if(envChange){
		environmentChange(Simulator::Now().GetSeconds());
	}
	m_statusTimer = Simulator::Schedule (Seconds(NEIGHBOR_CHECK_INTERVAL+interval), &MappleOnlineApplication::NeighborStatusTimer, this);

}
//call from sendSummaryToNeighbors()..
int MappleOnlineApplication::getNeighborRSSI(uint32_t id)
{
	return Neighbors[id].last_rssi;

}
//call from sendSummaryToNeighbors()..
int MappleOnlineApplication::getNeighborSNR(uint32_t id)
{
	return Neighbors[id].last_snr;

}
//1. call from printSamples()
uint8_t MappleOnlineApplication::getSendersNeighborsRssiRates(int32_t rssi){
	if(rssi > RSSI_RATE_LIMIT_MED){
		return RSSI_RATE_HIGH;
	}else if(rssi <= RSSI_RATE_LIMIT_MED && rssi > RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_MED;
	}else if(rssi < RSSI_RATE_LIMIT_LOW){
		return RSSI_RATE_LOW;
	}
	return 0;
}

/*
 * Called at time specified by Stop in a simulation script
 */
void MappleOnlineApplication::StopApplication ()
{
	NS_LOG_FUNCTION (this);
	printf("ApplicationStop at:%u now\n",myid);


	if (sendingSocket) {
		sendingSocket->Close ();
	}
	if(m_statusTimer.IsRunning()){
		m_statusTimer.Cancel();
		//  Simulator::DoScheduleDestroy(m_statustimer);
	}
	if(summaryTimer.IsRunning()){
		summaryTimer.Cancel();
		//	  Simulator::DoScheduleDestroy(summaryTimer);
	}
	if(event_estimator.IsRunning()){
		event_estimator.Cancel();
		//	  Simulator::DoScheduleDestroy(event_estimator);
	}
	//free the resources which we allocated in startapplication

//	fclose(lwpr_file);
	lwpr_free_model(&lwpr_model);
	//	free(shared_samples->sharedfeatures);
}
//I use three directive to get control over the code..
//a). NORMALIZEDSAMPLES (inorder to get normalized samples)
//b). SHARED_SAMPLES (which i use inside the normalized directive inorder to share the normalized features)
//c). RAWSAMPLES (Only for get raw samples without normalization)

//1. call from closeSample()
void MappleOnlineApplication::printSample(uint32_t id){

  //printf("\n-->PrintSamples at:%u at time:%lf for node:%u\n" ,GetNode()->GetId(),Simulator::Now().GetSeconds(),id);
	uint32_t i,sindex,sjix;
	int nodecount[10];
	int counter=1;
	double rssi=0.0;
	double snr=0.0;
	double rate=0;
	double prr=0.0;

#ifdef NORMALIZEDSAMPLES
	double avgrssi=0.0;
	double avgsnr = 0.0;
	double nodenormalized = 0.0;
	double snormalized= 0.0;
	double mcrate = 0.0;

#endif

	n_samples +=1;
	uint8_t rssi_rate = 0;
	int ncount=1;
	memset(&nodecount,0,sizeof(nodecount));
	//now we storing the normalized feature values
#ifdef RAWSAMPLES
	rssi = (1.0*(Samples[id].closed_sum_rssi) / (Samples[id].packets_recv_closed));
	writefeatures[1] = rssi;
	//	printf("The average RSSI is:%f\n",rssi);
#endif

#ifdef NORMALIZEDSAMPLES
	avgrssi = (1.0*Samples[id].closed_sum_rssi) / Samples[id].packets_recv_closed; //this is average value of the rssi for 50 packets and we needed this
	rssi = (avgrssi - sample_featuremin[1])/(sample_featuremax[1] - sample_featuremin[1]);//this is just last rssi value when we received the packet from this particular node
	writefeatures[1] = rssi;
	samples.samples_array[1] = rssi;
#endif//end of NORMALIZEDSAMPLES
	//	rssi = Samples[id].closed_sum_rssi / Samples[id].received; //this is average value of the rssi for 50 packets and we needed this
	//	snr = Samples[id].closed_sum_snr / Samples[id].received; // same for snr
#ifdef RAWSAMPLES
	snr = (1.0*Samples[id].closed_sum_snr) / Samples[id].packets_recv_closed; // same for snr
	writefeatures[2] = snr;
#endif

#ifdef NORMALIZEDSAMPLES
	avgsnr = (1.0*Samples[id].closed_sum_snr) / Samples[id].packets_recv_closed; // same for snr
	snr = (avgsnr - sample_featuremin[2]) /  (sample_featuremax[2] - sample_featuremin[2]);
	samples.samples_array[2] = snr;
	writefeatures[2] = snr;
#endif//end of NORMALIZEDSAMPLES

#ifdef RAWSAMPLES
	rate = Samples[id].sender_rate;
	writefeatures[3] = rate;
#endif

#ifdef NORMALIZEDSAMPLES
	rate = (Samples[id].sender_rate - sample_featuremin[3]) / (sample_featuremax[3] - sample_featuremin[3]); //a-->b , sending rate of node "a"
	samples.samples_array[3] = rate;
	writefeatures[3] = rate;
#endif//end of NORMALIZEDSAMPLES
	//prr is already normalized so we do not need to do it again
	prr = (1.0 * Samples[id].packets_recv_closed )/ (Samples[id].packets_sent_closed + 1);
	writefeatures[0] = prr;
	samples.samples_array[0] = prr;
	//	printf("Now:%f----------------->other values of prr is:%u\n",prr,Samples[id].packets_sent_closed/Samples[id].packets_recv_closed);
	//	printf("Packet_Sent_Closed:%u and Recv_Closed:%u\n",Samples[id].packets_sent_closed,Samples[id].packets_recv_closed);
	//#ifdef RAWSAMPLES
	//	fprintf(fp, "%6.2f 1:%f 2:%f 3:%f ", prr,rssi, snr, rate);
	//	printf("Printing raw samples\n");
	//#endif
	//check on given node id,
	for(i = 0; i < MAX_NODE_ID; i++) //MAX_NODE_ID //however I am Encountering the garbage values while sending 600
	{
		if( i == id)
			continue;
		if(is_set(Samples[id].receiver_status[i/32],i%32)){
			rssi_rate = getCurrentRssiRate(i);
			//snr = getNeighborSNR(i); // at the moment I am doing nothing with SNR rates of Receiver's Neighbors
			rate = Neighbors[i].neighbors_rate;
			for(uint32_t i=1;i<=3;i++){
				for(uint32_t j=1;j<=3;j++){
					if(rssi_rate == i && rate == j){
						nodecount[counter] = nodecount[counter] + 1;
						//go to the rssi + rate combination index and increase the number of node who have this particular combination
						//this is one solution, other is to write all 9 combination separatly in if statements/switch and  go.
					}
					counter++;
				}
			}
			//our own srrounding neighbors and their rates and rssi infos for this particular neighbor "id"!
			counter = 1;
		}
	}
	//there are total 9 combinations of RSSI[low,med,high] + RATE[low,med,high], and index of the nodecount presents these combination,
	//for example RSSI(low)+Rate(low) = 1, Rssi(low)+Rate(med) = 2, Rssi(low)+Rate(high) =3, Rssi(med)+Rate(low) and so on..!
	//and contents of the array hold the Number of neighbors how have these rate combinations.!
	for(ncount =1; ncount<=9 ; ncount++){
#ifdef NORMALIZEDSAMPLES
		nodenormalized = (nodecount[ncount] - sample_featuremin[ncount+3]) / (sample_featuremax[ncount+3] - sample_featuremin[ncount+3]);
		writefeatures[ncount+3] = nodenormalized;
		//fprintf(fp,"%d:%f ",ncount+3,nodenormalized);
		samples.samples_array[ncount+3] = nodenormalized;
#endif//end of NORMALIZEDSAMPLES

#ifdef RAWSAMPLES
		writefeatures[ncount+3] = nodecount[ncount];
		//		fprintf(fp,"%d:%d ",ncount+3,nodecount[ncount]);
#endif
	} // end of Receiver Status loop

	int snodecount[10];
	int scounter = 1;
	int scount = 1;
	memset(&snodecount,0,sizeof(snodecount));
	//fprintf(fp,"] sender_status [");
	sjix=0;
	for(sindex = 0; sindex < MAX_NODE_ID; sindex++)
	{
		if(is_set(Samples[id].sender_status[sindex/32],sindex%32)){
			rssi = Samples[id].sender_rssi_data[sjix];
			//	snr = Samples[id].sender_snr_data[sjix]; // doins nothing with Sender's neighbors Snr rates.. see them later!
			rate = Samples[id].sender_rate_data[sjix];
			for(uint32_t i=1;i<=3;i++){
				for(uint32_t j=1;j<=3;j++){

					if(getSendersNeighborsRssiRates(rssi) == i && rate == j){
						snodecount[scounter] = snodecount[scounter] + 1;
						//go to the rssi + rate combination index and increase the number of node who have this particular combination
						//this is one solution, other is to write all 9 combination separatly in if statements/switch and  go.
					}
					scounter++;
				}
			}
			scounter = 1;
			sjix=sjix+1;
		}
	}

	for(scount =1; scount<=9 ; scount++){

#ifdef NORMALIZEDSAMPLES
		snormalized = (snodecount[scount] - sample_featuremin[scount+12]) / (sample_featuremax[scount+12] - sample_featuremin[scount+12]); //sample_feature[] starting from zero so
		writefeatures[scount+12] = snormalized;
		samples.samples_array[scount+12] = snormalized;
		//fprintf(fp,"%d:%f ",scount+12,snormalized);
#endif //end of NORMALIZEDSAMPLES

#ifdef RAWSAMPLES
		//		fprintf(fp,"%d:%d ",scount+12,snodecount[scount]);
		writefeatures[scount+12] = snodecount[scount];
#endif
		//TSN (total sender's Neighbors nodes who holding this particular rates combinations).. see the description in above
		//receiver scenarios.!
	}  /// end of Sender Status loop

#ifdef NORMALIZEDSAMPLES
	mcrate = (getCurrentSendingRate() - sample_featuremin[22]) / (sample_featuremax[22] - sample_featuremin[22]);
	writefeatures[22] = mcrate;
	samples.samples_array[22] = mcrate;

	double distance = (getSenderDistance(id) - sample_featuremin[23]) / (sample_featuremax[23] - sample_featuremin[23]);
	writefeatures[23] = distance;
	samples.samples_array[23] = distance;

//	writefeatures[23] = getSenderDistance(id);
//	samples.samples_array[23] = getSenderDistance(id); // at the moment this unormalized values

	//after last sample we will write into node sample files..!

	//for sake of safty we are storing the samples and test feature inside a queue inorder to tolerate the faults
	// later the separate timer will hande these queue to pop the elements. During samples collection we may
	//suffer access of samples due to high traffice to cop this loss we are using queue.(Actually I am wishing for
	//deque where they give us the more control over the data.
	samples_queue.push(samples); // build a queue which hold all incomming samples
	writeSamplesFile();
//	printf("-------------------Sample Q size:%d \n",samples_queue.size());
	//fprintf(fp,"22:%f\n",mcrate);
#endif  //end of normalized directive




#ifdef RAWSAMPLES
	//	fprintf(fp,"22:%u\n",getCurrentSendingRate());
	writefeatures[22] = getCurrentSendingRate();
	writefeatures[23] = getSenderDistance(id);
	memcpy(samples.samples_array,writefeatures,sizeof(double)*24);
	samples_queue.push(samples);
	writeSamplesFile();

#endif
	//for sharing the samples with other nodes we need to invoke this directive.!
#ifdef SHARED_SAMPLES
	//memcpy(recentSamples,writefeatures,sizeof(writefeatures));
	//here we need separate schedualar who will send the samples to other nodes
	sharedSamples(writefeatures);
	//ssampleEvent =Simulator::ScheduleNow(&MappleOnlineApplication::sharedSamples,this);
#endif


} // End of printSample() function



//1. call from the printSamples()
double MappleOnlineApplication::getSenderDistance(uint32_t senderid){


	Vector myposition = m_mobility->GetPosition();

	Ptr<Node> senderNode = m_nodeContainer.Get(senderid);

	Ptr<RandomWaypointMobilityModel> sendermobility = senderNode->GetObject<RandomWaypointMobilityModel>();

	Vector senderposition = sendermobility->GetPosition();
	Vector difference;

	difference.x = myposition.x - senderposition.x;
	difference.y = myposition.y - senderposition.y;
	difference.z = myposition.z - senderposition.z;
	//here we just find the Euclidean distance from the sender node to the node who receive this sample
	double e_distance = sqrt(pow(difference.x,2)+pow(difference.y,2)+pow(difference.z,2));

	//printf("Node %u far from node %u[x:%f,y:%f,z:%f] edistance:%f\n",myid,senderid,difference.x,difference.y,difference.z,e_distance);

	return e_distance;

}


///start of sample comparisions between two vectors
#ifdef SHARED_SAMPLES
//1. Call from printSamples()
//instead of passing all shared_Sample structure I will only pass the sample values which want to compare?
void MappleOnlineApplication::sharedSamples(double recentsamples[]){

	//	printf("Now in sharedSamples fucn: %f %f %f %f \n",recentsamples[0],recentsamples[1],recentsamples[2],recentsamples[3]);
	double type = 0;
	//before we procceed to send the samples we will make little assessment for the quality of the sample.
	//we will check either we sending the same sample which we received as previoused or new one
	if(n_samples == 1){
		//		oldSamples = recent;
		//	printf("sharedSamples n_samples:%u\n",n_samples);
		memcpy(oldSamples,recentsamples,sizeof(double)*24);
	}
	if(n_samples > 1){
		type = SampleQualityTest(recentsamples,oldSamples);
		if(type < 1){
			//If samples are uncorelated OR negatively corelated.. send them
		//	printf("Sample is uncorelated :%f\n",type);
			sendSharedSamples(recentsamples);
		}
		memcpy(oldSamples,recentsamples,sizeof(double)*24);//after checking the 2nd sample we need it for third one to copmare it
	}
}
//Correlation r= N*Sum(xy)-(Sum(x))(Sum(y)) / sqrt[N*Sum(x^2) - (Sum(x))^2][N*Sum(y^2)-(Sum(y))^2]
//1. call from sharedSamples()
double MappleOnlineApplication::SampleQualityTest(double recentVector[],double oldVector[]){

	//find the corealtion formulat and find the difference between two sample vecotrs
	double r = 0;
	int N=24;
	double Sxy = 0;
	double Sx = 0;
	double Sy = 0;
	double Sx2 = 0;
	double x2 = 0;
	double y2 = 0;
	double Sy2 = 0;
	for(uint32_t i=0;i<=23;i++){

		Sx = Sx + recentVector[i];
		Sy = Sy + oldVector[i];

		x2 = pow(recentVector[i],2);
		y2 = pow(oldVector[i],2);

		Sx2 = Sx2 + x2;
		Sy2 = Sy2 + y2;

		Sxy = Sxy + (recentVector[i]*oldVector[i]);

	}

	r = (N*(Sxy) - ((Sx)*(Sy))) / sqrt((N*(Sx2) - pow(Sx,2))*(N*(Sy2) - pow(Sy,2)));
	return r;
}
//1 call from sharedSamples
void MappleOnlineApplication::sendSharedSamples(double Samples[]){

	shared_samples_t shared_samples;
	shared_samples.msg_type = SHAREDSAMPLES;
	shared_samples.node_id = myid;
	//loop to fill the svm_node with share_saple_t structure to send
	//	for(uint32_t count=0;count <= 22;count++){
	//
	//		shared_samples.sharedfeatures[count]= recentSamples[count];
	//	}
	//	printf("Sizeof share_Samples_t:%d\n",sizeof(shared_samples_t));

	memcpy(shared_samples.sharedfeatures,Samples,sizeof(double)*24);
	Ptr<Packet> packet = Create<Packet> (reinterpret_cast<const uint8_t*>(&shared_samples),sizeof(shared_samples_t));
	//destination = Ipv4Address::GetBroadcast();
	//	printf("Packet send Size is-------------%d\n",packet->GetSize());
	sharedSocket -> SendTo(packet,0,InetSocketAddress(sharedDestination,8000)); //
	//printf("Broadcast this sample\n");

}

#endif // end of shared sample functionality block

//writing in the file make no difference either the sample is as Shared sample or normal sample..
//inorder to use these samples for model building we are bit independent to how fast we receive, so
//in general this is would required seprate timer to deal with...

//but here is one trade off. that is we have samples in queue, we first write the sample and then satart training

//1. call from printSamples (the main function)
//2. call from the readSahredSamples()
bool MappleOnlineApplication::writeSamplesFile(){

	//	printf("we writing file now\n");

	if(!mutex){

		mutex = true;
//		double *q_data;// = (double*) malloc((23)*sizeof(double));
		samples_t writesamples;
	//	double q_data[23];
#ifdef TRAINING_PREDICTION
		double target_label=0.0;
#endif

		int size = samples_queue.size();

		while(size != 0){

		//	memcpy(q_data,samples_queue.front(),sizeof(double)*23);
			//	printf("Sizeof the qu in writesamplefiles is:%d\n",size);
			writesamples = samples_queue.front();//access items front of the queue..and write them in the file
			for(uint32_t count=0; count<=23;count++){

				if(count == 0){
//							printf("Write:%6.2f",writesamples.samples_array[0]);
#ifdef WRITE_SAMPLES
					fprintf(fp,"%6.2f ",writesamples.samples_array[0]);
#endif

#ifdef LWPR_TP
				//	printf("hereis problem\n");
			lwpr_y = writesamples.samples_array[0];
#endif

#ifdef TRAINING_PREDICTION
					target_label = writesamples.samples_array[0];// target lable to calculate the error..
#endif

				}else{
//						printf(" %u:%f",count,writesamples.samples_array[count]);
#ifdef WRITE_SAMPLES
					fprintf(fp,"%u:%f ",count,writesamples.samples_array[count]);
#endif

#ifdef LWPR_TP
			lwpr_x[count-1]  = writesamples.samples_array[count];
#endif

#ifdef TRAINING_PREDICTION
					testfeatures[count-1].index = count; //testfeatures index for livsvm formate
					testfeatures[count-1].value = writesamples.samples_array[count]; //testfeatures values to predict

//					printf("Testfeature values:%f\n",testfeatures[count-1]);
#endif
				}
			}
//					printf("\n");
#ifdef WRITE_SAMPLES
			fprintf(fp,"\n");
			//memset(q_data,0x0,sizeof(double)*23);
#endif

			//printf("its should be here\n");
#ifdef LWPR_TP
			//printf("Its now in lwpr tp:\n");
			Lwpr_TrainPredict(lwpr_x,lwpr_y);
#endif
			/////////////////////////////////////////////Here the sample training part Started///////////////////////////////////////
			//after writting file, we pass that file name into our Svm_Training() function with model name of this particular node
			//for training functionality I am using the seprate directive, might be at some time we don't need to using training
#ifdef TRAINING_PREDICTION
			if(n_samples == 1){
				Svm_Training(saved_sample.c_str(),model_name.c_str());
			}
			if(n_samples > 1 ){
				if((model=svm_load_model(read_model.c_str()))==0)
				{
					printf("[ERROR] can't open at sample # %d SVM model file: %s\n",n_samples,read_model.c_str());
					exit(1);
				}else{
					//			nooffeatures = model-> vecDim;
					//			//nooffeatures = 22;
					//			svm_type=svm_get_svm_type(model);
					//			nr_class=svm_get_nr_class(model);

					double predict_label = svm_predict(model,testfeatures);
					double error = abs(target_label - predict_label);
					if(!isnan(error)){
						errSum += (error * error); //Square the error and commulate all those values
						errCounter++; //how many data points collected so far
						fprintf(fp_error,"%u	%g\n",myid,errSum/errCounter);
						//printf("Error: %g ErrSum:%g errCounter:%u - predicted label:%g cumMSE:%g \n",error,errSum,errCounter,predict_label,errSum/errCounter);
					}
					Svm_Training(saved_sample.c_str(),model_name.c_str());
					//now again training on the same samples
					//here we need some schedualor for writeSamplesFile() in loop?
				}
				//Copy the link quality features into 1D array for testing with the SVM
				//Here I just alocate the memory for one time (no need to allocate every time), and next time we just overwrite the previous sample
				//		testfeatures = (struct svm_node*) malloc((nooffeatures+1)*sizeof(struct svm_node));
			}

#endif //end of TRAINING_PREDICTION directive
			samples_queue.pop(); // whatever we wrote we just delete that sample
			size = samples_queue.size();
		}//end of first if(size !=0 )

		mutex = false;

	}
	return true;
	//	Simulator::ScheduleNow(&MappleOnlineApplication::writeSamplesFile, this);
}
///////////////////////////////////////LWPR TRAINING AND PREDICTION BEGIN///////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef LWPR_TP
void MappleOnlineApplication::Lwpr_TrainPredict(double datax[], double datay){


	//printf("yes lpwr_trainpredict called\n");
	double ypredict;
	double ytrain;
	if(lwpr_counter == 0){
		lwpr_update(&lwpr_model, datax, &datay, &ytrain, NULL);
	lwpr_mse_T +=  (datay-ytrain)*(datay-ytrain);
	lwpr_counter++;
	}else{

		lwpr_predict(&lwpr_model, datax, 0.001, &ypredict, NULL, NULL);
//		ypredict = lwpr_model.predict(datax);
		lwpr_mse_P += (datay-ypredict)*(datay-ypredict);
		lwpr_update(&lwpr_model, datax, &datay, &ytrain, NULL);;
		lwpr_mse_T +=  (datay-ytrain)*(datay-ytrain);
		lwpr_counter++;
		//printf("Writing lwpr errors\n");
		fprintf(lwpr_file,"%u	%f  %f\n",myid,lwpr_mse_T/lwpr_counter, lwpr_mse_P/lwpr_counter);


	}


}

#endif
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

///////////////////////////////////////SVM-TRAIN.C BEGING///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef TRAINING_PREDICTION
void print_null(const char *s) {}

void exit_with_help()
{
	printf(
			"Usage: svm-train [options] training_set_file [model_file]\n"
			"options:\n"
			"-s svm_type : set type of SVM (default 0)\n"
			"	0 -- C-SVC\n"
			"	1 -- nu-SVC\n"
			"	2 -- one-class SVM\n"
			"	3 -- epsilon-SVR\n"
			"	4 -- nu-SVR\n"
			"-t kernel_type : set type of kernel function (default 2)\n"
			"	0 -- linear: u'*v\n"
			"	1 -- polynomial: (gamma*u'*v + coef0)^degree\n"
			"	2 -- radial basis function: exp(-gamma*|u-v|^2)\n"
			"	3 -- sigmoid: tanh(gamma*u'*v + coef0)\n"
			"	4 -- precomputed kernel (kernel values in training_set_file)\n"
			"-d degree : set degree in kernel function (default 3)\n"
			"-g gamma : set gamma in kernel function (default 1/num_features)\n"
			"-r coef0 : set coef0 in kernel function (default 0)\n"
			"-c cost : set the parameter C of C-SVC, epsilon-SVR, and nu-SVR (default 1)\n"
			"-n nu : set the parameter nu of nu-SVC, one-class SVM, and nu-SVR (default 0.5)\n"
			"-p epsilon : set the epsilon in loss function of epsilon-SVR (default 0.1)\n"
			"-m cachesize : set cache memory size in MB (default 100)\n"
			"-e epsilon : set tolerance of termination criterion (default 0.001)\n"
			"-h shrinking : whether to use the shrinking heuristics, 0 or 1 (default 1)\n"
			"-b probability_estimates : whether to train a SVC or SVR model for probability estimates, 0 or 1 (default 0)\n"
			"-wi weight : set the parameter C of class i to weight*C, for C-SVC (default 1)\n"
			"-v n: n-fold cross validation mode\n"
			"-q : quiet mode (no outputs)\n"
	);
	exit(1);
}

void exit_input_error(int line_num, char check)
{
	fprintf(stderr,"Wrong input format at line %d:%c\n", line_num,check);
	exit(1);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name);
void read_problem(const char *filename);
void do_cross_validation();

struct svm_parameter param;		// set by parse_command_line
struct svm_problem prob;		// set by read_problem
struct svm_model *model;
struct svm_node *x_space;
int cross_validation;
int nr_fold;

static char *line = NULL;
static int max_line_len;

static char* readline1(FILE *input)
{
	int len;

	if(fgets(line,max_line_len,input) == NULL)
		return NULL;

	while(strrchr(line,'\n') == NULL)
	{
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if(fgets(line+len,max_line_len-len,input) == NULL)
			break;
	}
	return line;
}
// we change the main fucntion to Svm_Training
//int main(int argc, char **argv)
//1. call from Print_Samples();
int MappleOnlineApplication::Svm_Training(const char *inputfile,const char *modelfile)
{
	char input_file_name[1024];
	char model_file_name[1024];
	const char *error_msg;
	strcpy(input_file_name,inputfile);
	strcpy(model_file_name,modelfile);
	int argc = 0;
	char **argv = NULL;
	//	printf("input file name is: %s\n",input_file_name);
	//	printf("output modelname is: %s\n",model_file_name);

	parse_command_line(argc, argv, input_file_name, model_file_name);
	read_problem(input_file_name);
	error_msg = svm_check_parameter(&prob,&param);

	if(error_msg)
	{
		fprintf(stderr,"ERROR: %s\n",error_msg);
		exit(1);
	}

	if(cross_validation)
	{
		do_cross_validation();
	}
	else
	{
		model = svm_train(&prob,&param);
		if(svm_save_model(model_file_name,model))
		{
			fprintf(stderr, "can't save model to file %s\n", model_file_name);
			exit(1);
		}
		svm_free_and_destroy_model(&model);
	}
	svm_destroy_param(&param);
	free(prob.y);
	free(prob.x);
	free(x_space);
	free(line);

	return 0;
}

void do_cross_validation()
{
	int i;
	int total_correct = 0;
	double total_error = 0;
	double sumv = 0, sumy = 0, sumvv = 0, sumyy = 0, sumvy = 0;
	double *target = Malloc(double,prob.l);

	svm_cross_validation(&prob,&param,nr_fold,target);
	if(param.svm_type == EPSILON_SVR ||
			param.svm_type == NU_SVR)
	{
		for(i=0;i<prob.l;i++)
		{
			double y = prob.y[i];
			double v = target[i];
			total_error += (v-y)*(v-y);
			sumv += v;
			sumy += y;
			sumvv += v*v;
			sumyy += y*y;
			sumvy += v*y;
		}
		printf("Cross Validation Mean squared error = %g\n",total_error/prob.l);
		printf("Cross Validation Squared correlation coefficient = %g\n",
				((prob.l*sumvy-sumv*sumy)*(prob.l*sumvy-sumv*sumy))/
				((prob.l*sumvv-sumv*sumv)*(prob.l*sumyy-sumy*sumy))
		);
	}
	else
	{
		for(i=0;i<prob.l;i++)
			if(target[i] == prob.y[i])
				++total_correct;
		printf("Cross Validation Accuracy = %g%%\n",100.0*total_correct/prob.l);
	}
	free(target);
}

void parse_command_line(int argc, char **argv, char *input_file_name, char *model_file_name)
{
	int i;
	void (*print_func)(const char*) = NULL;	// default printing to stdout

	// default values
	param.svm_type = EPSILON_SVR;
	param.kernel_type = RBF;
	param.degree = 3;
	param.gamma = 0.01;	// 1/num_features
	param.coef0 = 0;
	param.nu = 0.5;
	param.cache_size = 100;
	param.C = 1;
	param.eps = 1e-3;
	param.p = 0.1;
	param.shrinking = 1;
	param.probability = 0;
	param.nr_weight = 0;
	param.weight_label = NULL;
	param.weight = NULL;
	cross_validation = 0;

	// parse options
	for(i=1;i<argc;i++)
	{
		if(argv[i][0] != '-') break;
		if(++i>=argc)
			exit_with_help();
		switch(argv[i-1][1])
		{
		case 's':
			param.svm_type = atoi(argv[i]);
			break;
		case 't':
			param.kernel_type = atoi(argv[i]);
			break;
		case 'd':
			param.degree = atoi(argv[i]);
			break;
		case 'g':
			param.gamma = atof(argv[i]);
			break;
		case 'r':
			param.coef0 = atof(argv[i]);
			break;
		case 'n':
			param.nu = atof(argv[i]);
			break;
		case 'm':
			param.cache_size = atof(argv[i]);
			break;
		case 'c':
			param.C = atof(argv[i]);
			break;
		case 'e':
			param.eps = atof(argv[i]);
			break;
		case 'p':
			param.p = atof(argv[i]);
			break;
		case 'h':
			param.shrinking = atoi(argv[i]);
			break;
		case 'b':
			param.probability = atoi(argv[i]);
			break;
		case 'q':
			print_func = &print_null;
			i--;
			break;
		case 'v':
			cross_validation = 1;
			nr_fold = atoi(argv[i]);
			if(nr_fold < 2)
			{
				fprintf(stderr,"n-fold cross validation: n must >= 2\n");
				exit_with_help();
			}
			break;
		case 'w':
			++param.nr_weight;
			param.weight_label = (int *)realloc(param.weight_label,sizeof(int)*param.nr_weight);
			param.weight = (double *)realloc(param.weight,sizeof(double)*param.nr_weight);
			param.weight_label[param.nr_weight-1] = atoi(&argv[i-1][2]);
			param.weight[param.nr_weight-1] = atof(argv[i]);
			break;
		default:
			fprintf(stderr,"Unknown option: -%c\n", argv[i-1][1]);
			exit_with_help();
		}
	}

	svm_set_print_string_function(print_func);

	// determine filenames

	//	if(i>=argc)
	//		exit_with_help();
	//
	//	strcpy(input_file_name, argv[i]);
	//
	//	if(i<argc-1)
	//		strcpy(model_file_name,argv[i+1]);
	//	else
	//	{
	//		char *p = strrchr(argv[i],'/');
	//		if(p==NULL)
	//			p = argv[i];
	//		else
	//			++p;
	//		sprintf(model_file_name,"%s.model",p);
	//	}
}

// read in a problem (in svmlight format)

void read_problem(const char *filename)
{
	int elements, max_index, inst_max_index, i, j;
	FILE *fp = fopen(filename,"r");
	char *endptr;
	char *idx, *val, *label;

	if(fp == NULL)
	{
		fprintf(stderr,"can't open input file %s\n",filename);
		exit(1);
	}

	prob.l = 0;
	elements = 0;

	max_line_len = 1024;
	line = Malloc(char,max_line_len);
	while(readline1(fp)!=NULL)
	{
		char *p = strtok(line," \t"); // label

		// features
		while(1)
		{
			p = strtok(NULL," \t");
			if(p == NULL || *p == '\n') // check '\n' as ' ' may be after the last feature
				break;
			++elements;
		}
		++elements;
		++prob.l;
	}
	rewind(fp);

	prob.y = Malloc(double,prob.l);
	prob.x = Malloc(struct svm_node *,prob.l);
	x_space = Malloc(struct svm_node,elements);

	max_index = 0;
	j=0;
	for(i=0;i<prob.l;i++)
	{
		inst_max_index = -1; // strtol gives 0 if wrong format, and precomputed kernel has <index> start from 0
		readline1(fp);
		prob.x[i] = &x_space[j];
		label = strtok(line," \t\n");
		if(label == NULL) // empty line
			exit_input_error(i+1,'a');

		prob.y[i] = strtod(label,&endptr);
		if(endptr == label || *endptr != '\0')
			exit_input_error(i+1,'b');


		while(1)
		{
			idx = strtok(NULL,":");
			val = strtok(NULL," \t");

			if(val == NULL)
				break;

			errno = 0;
			x_space[j].index = (int) strtol(idx,&endptr,10);
			if(endptr == idx || errno != 0 || *endptr != '\0' || x_space[j].index <= inst_max_index)
				exit_input_error(i+1,'c');
			else
				inst_max_index = x_space[j].index;

			errno = 0;
			x_space[j].value = strtod(val,&endptr);
			if(endptr == val || errno != 0 || (*endptr != '\0' && !isspace(*endptr)))
				exit_input_error(i+1,'d');

			++j;
		}

		if(inst_max_index > max_index)
			max_index = inst_max_index;
		x_space[j++].index = -1;
	}

	if(param.gamma == 0 && max_index > 0)
		param.gamma = 1.0/max_index;

	if(param.kernel_type == PRECOMPUTED)
		for(i=0;i<prob.l;i++)
		{
			if (prob.x[i][0].index != 0)
			{
				fprintf(stderr,"Wrong input format: first column must be 0:sample_serial_number\n");
				exit(1);
			}
			if ((int)prob.x[i][0].value <= 0 || (int)prob.x[i][0].value > max_index)
			{
				fprintf(stderr,"Wrong input format: sample_serial_number out of range\n");
				exit(1);
			}
		}

	fclose(fp);
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
/////////////////////////////////////////////////End of NS3 NameSpace boundary/////////////////////////////////////////////
} // end of namespace ns3
