/*
 * MAPPLE-ONLINE-aap.h
 *
 *  Created on: Jul 16, 2012
 *      Author: imran
 */

#ifndef MAPPLE_ONLINE_AAP_H_
#define MAPPLE_ONLINE_AAP_H_

#include "ns3/application.h"
#include "ns3/event-id.h"
#include "ns3/ptr.h"
#include "ns3/traced-callback.h"
#include "ns3/address.h"
#include "ns3/ipv4-interface-container.h"
#include "ns3/node-container.h"
#include "ns3/net-device.h"
#include "ns3/random-variable.h"
#include "ns3/mobility-model.h"
#include "ns3/waypoint-mobility-model.h"
#include "ns3/random-waypoint-mobility-model.h"
#include "ns3/waypoint.h"
#include "ns3/wifi-mac-queue.h"
#include <iostream>
#include <string>
#include <map>
#include <queue>

#include <stdio.h>

#include <cmath>

#include "svm.h" // header for SVM library

//#include "lwpr_cc.h" // header for lwpr library

#include "lwpr.h"



namespace ns3{

class Socket;
class Packet;
class Address;

//#define LWPR_TP 1  			//inorder to use LWPR regressor uncoment this directive
#define SHARED_SAMPLES 1   //inorder to share sample for cooperative learning please uncoment this directive
//#define RAWSAMPLES 1    //printing raw samples without normalization please uncomment this
#define WRITE_SAMPLES 1 //while using lwpr stuff you need to comment this directive, otherwise keep uncomment for RAWSAMPLES and SVM STUFF...!
#define NORMALIZEDSAMPLES 1 //make our samples normalized
#define TRAINING_PREDICTION 1 //for turning SVM training and prediction for this we also need above directive to turn on

#define SHAREDSAMPLES 23

#define START_TRAINING 0
#define RATE_ESTIMATION_INTERVAL 4096 //instead giving it static value I need to put it in typeid accessor where script decide the vale
#define LQE_RATE_LIMIT_MED 1650 //orginally Eduardo rate is 1500
#define LQE_RATE_LIMIT_HIGH  3300 //orignally 3000

#define MIN_PACKETS_CLOSED 20
#define MAX_PACKETS_CLOSED 150

#define LQE_RATE_LOW  1
#define	LQE_RATE_MED  2
#define	LQE_RATE_HIGH 3

#define RSSI_RATE_LIMIT_MED  -50 // IF THE Rssi value is greater than this we say its medium otherwise High
#define RSSI_RATE_LIMIT_LOW  -80 // Greater than this, mean Low Rssi
//At the moment we just create 3 bins(crossponding to 3 features), //infuture we might increase
#define RSSI_RATE_LOW  1
#define	RSSI_RATE_MED  2
#define	RSSI_RATE_HIGH 3


#define MAX_NODE_ID 100
#define MAX_REPORTED_NEIGHBORS 50
#define SEND_UNICAST 0 //I really dont know why Michal use this macro?
#define LQE_RATE_LOW 1
#define	LQE_RATE_MED 2
#define	LQE_RATE_HIGH 3
#define SUMMARY_PACKET 9

#define STATE_NOT_ACTIVE 0
#define	STATE_ACTIVE  1

#define IMMEDIATE_SUMMARY 0

#define MAX_NEIGHBOR_TTL 3.0  // seconds
#define NEIGHBOR_CHECK_INTERVAL 2.048 // seconds //2048 milliseconds


// precomputed normalization values (SVM PREDICTION)
static const double sample_featuremax[] = {0.993631,18,111,3,7,6,7,4,5,5,2,1,2,7,6,7,4,5,5,2,1,2,3,126.7609};
static const double sample_featuremin[] = {0.040000,-91,2,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1.3303};



////////////////these structures copied from the Eduardo code for online Application////////////////////
//we different names otherwise it would conflict with MAPPLE-app.h definitions ..
typedef struct onlinesample {
	// Flag indicating is sample is "open"
	bool open;
	// Bit array indicating active neighbors at sender nod
	uint32_t sender_status[5];
	uint32_t receiver_status[5];

	int8_t sender_rssi_data[MAX_REPORTED_NEIGHBORS];
	uint8_t sender_rate_data[MAX_REPORTED_NEIGHBORS]; // this senders neighbors data rates
	uint8_t sender_snr_data[MAX_REPORTED_NEIGHBORS]; //At the moment we also consider this.

	uint8_t sender_rate;   //sender's own applications sending rates, however source don't need to consider its sending rates at the moment
	uint8_t activeNeighbors; // number of active neighbors of the source

	uint32_t first_counter; // just need to track the first summary packet, which we filled at openSample() and it remains that, we will not update it, unless we received new sampel
	uint32_t packets_recv_closed; //using below received value we count how many packet received
	uint32_t packets_sent_closed; //for the receiver node, sake of statistics we count how many packets we sent

	uint32_t sent_traffic_counter_lastclose; //
	int32_t last_rssi;
	uint32_t received;  //just to indicated that we have received one packet, kind a fix value al the time
	int32_t closed_sum_rssi;
	int32_t closed_sum_snr;

	int32_t sum_rssi; // this is just for calculating the average for the particular node which we have connection at the moment..a->b
	int32_t sum_snr;  //orignally its not there, but I am going to add it!

	double open_time;  //in openSample() we filled with current simulation time, remain same time for future updated same sample
	double last_close_time; //for one particular sample, when we received it last time, usually current local time stamp
	//Michal's Additional features which not supported by Eduardo implementation
	//at the moment not needed
	//	double last_sample_time;
	//	double avgTravelTime; // current_time - sent_time of the summary message might be helpful feature, but could be unperdictable
	//	double avgSendTime; //Michal! calculated starting from the moment when dca-txop requested for channel access for this packet...
	//	double avgSendTimeWithoutBackoff; //Michal! // calculated starting from the moment when packet was dequeued from the MAC queue..
	//	uint32_t queueSize; //MAC queue size is a important feature, to decide routing paths..
} online_sample_t;

typedef struct summary_msg{
	uint8_t  msg_type;
	//	uint16_t sender_id;
	uint32_t source_status[5];
	uint32_t last_update;  //Last time node state was updated
	uint8_t my_current_rate; //this is for sender owns overall traffic rate need for PRR.
	uint8_t activeNeighbors;
	int8_t neighbors_rssi_data[MAX_REPORTED_NEIGHBORS]; //this for neighbors rssi information
	uint8_t neighbors_snr_data[MAX_REPORTED_NEIGHBORS]; //this for neighbors rssi information
	uint8_t neighbors_rate_data[MAX_REPORTED_NEIGHBORS]; //this for neighbors sendding rates information..
	uint32_t counter;
	uint32_t sent_traffic_counter; //since last summary message counter
	//	double sendTime;//Michal Feature inorder to calculate the timings
}summary_msg_t;

typedef struct online_neighbor_entry {
	float distance;  //at the moment not calculating
	int32_t last_rssi; //same here
	int32_t last_snr; //last SNR ratio of the source node traffic, not nessarily summary packet
	uint8_t state; //state active or not
	uint8_t neighbors_rate;
	uint64_t last_update; //when we receive last time summary message, local time of the source node summary packet generation
	double last_message_time;
} online_neighbor_entry_t;

typedef struct online_rate_entry {
	int32_t bytes_sent;
	uint32_t rate;  // in bytes per second
	float packet_len;
	float interval;
} online_rate_entry_t;

//Inorder to share samples we need to define the particular packet for this particular need

typedef struct shared_samples {

	uint8_t msg_type;
	uint32_t node_id;
	//before it was 23 three features now we add the euclidean distance, and total features become 24
	double sharedfeatures[24];


}shared_samples_t;

typedef struct samples{
	double samples_array[24];
}samples_t;

////////////////////////End of Structures Definitions/////////////////////
#define is_set(X,Y) (X & (1 << (Y)))
#define set_bit(X,Y) X = ((X) | (1 << (Y)))
#define clear_bit(X,Y) X = ((X) & ~(1 << (Y)))

#define Malloc(type,n) (type *)malloc((n)*sizeof(type))
typedef std::vector<double> doubleVec;
//lwpr Macros ..
#define SEED_RAND()     srand48(time(NULL))
#define URAND()         drand48()

///////////////Start Definition of the MappleOnlineApplication class ..////////////////////////
class MappleOnlineApplication : public Application
{


public:

	static TypeId GetTypeId (void); //this is for typeId system of ns3, we must define this in-order to make application accessible
	MappleOnlineApplication ();   // constructor define our own
	virtual ~MappleOnlineApplication();
	uint32_t GetTotalRx (void) const; //get total bytes received

	Ptr<Socket> GetListeningSocket (void) const;

	/**
	 * \return list of pointers to accepted sockets
	 */
	std::list<Ptr<Socket> > GetAcceptedSockets (void) const;
	///inline function
	void SetInterfaceContainer(Ipv4InterfaceContainer inter){this->interfaces = inter;}
	//the helper call this function and set the node container..and from here we can access the id and addresses
	void SetNodeContainer(NodeContainer nc){this->m_nodeContainer = nc;}

	void SetSimulationtime(double simulationtime){this->m_Simulationtime = simulationtime; }

	void SetPrefix(std::string prefix) {this->outFilePath = prefix;}
	//static void SendingRate(void);
	void RateEstimationTimer(void);
	void sendSummaryToNeighbors();

	//void SetModel(std::string model) {this->modelName = model;}
	void environmentChange(double);
private:

	virtual void StartApplication (void);    // Called at time specified by Start
	virtual void StopApplication (void);    // Called at time specified by Start
	void readSharedSamples (Ptr<Socket>); //callback as handleread
	void HandleAccept (Ptr<Socket>, const Address& from);
	void HandlePeerClose (Ptr<Socket>);
	void HandlePeerError (Ptr<Socket>);

	bool HandlePromiscuous(Ptr<NetDevice>, Ptr<const Packet>, uint16_t, const Address &, const Address &, enum NetDevice::PacketType);	// Called for packets successfully received by the device, in promiscuous mode
	bool HandleSendTraffic(const Address &, uint32_t, double , uint16_t);	// Called for packets successfully received by the device, in promiscuous mode
	uint8_t getCurrentRssiRate(uint32_t); //neighbor[].last_rssi would help us to find the current rssi rate of sender node
	uint8_t getOldRssiRate(uint32_t);
	uint8_t getSendersNeighborsRssiRates(int32_t);
	uint8_t getCurrentSendingRate(void);
	bool isNeighborStatusActive(uint32_t);
	bool isFirstPacketSample(uint32_t);
	void BuildIdMapp(void); //this is mapp for Ip address to id's
	bool setNeighborStatus(uint32_t id, uint8_t state);
	bool myNeighborsRateChanged(summary_msg_t *, uint32_t id);



	void handlePacketReception(summary_msg_t*, uint32_t, double, double);
	void processNeighborSummary(summary_msg_t*, uint32_t, double, double);

	void openSamplesDatabase(summary_msg_t*, uint32_t, double, double);
	void updateSamplesDatabase(summary_msg_t*, uint32_t, double, double);
	void printSample(uint32_t id);

	bool hasOpenSample(uint32_t);
	void closeSample(uint32_t);
	void closeAllSamples(void);

	bool isNeighborEnvironmentChanged(summary_msg_t *, uint32_t);

	void NeighborStatusTimer(void);

	double getSenderDistance(uint32_t); // Get sender distance whome we are saving the sample

	int getNeighborRSSI(uint32_t);
	int getNeighborSNR(uint32_t);


	//function for sharing the samples
	bool writeSamplesFile();

	void sharedSamples(double []);
	double SampleQualityTest(double [],double []);
	void sendSharedSamples(double []);

	//struct svm_model * StartTraining(uint32_t id);
	int Svm_Training(const char *,const char *);

	//lwpr training and predictions.
	void Lwpr_TrainPredict(double [],double);

	//	void read_problem(const char *);
	//tracedcallback functions just used incontrast to collect some inner states of the application..
	TracedCallback<Ptr<const Packet>, const Address &> m_RecvPacketTrace;		// a trace source to detect received packages
	TracedCallback<Ptr<const Packet> > m_TransPacketTrace;						// a trace source to detect transmitted packages

	uint16_t portBase;
	TypeId  m_tid;          // Protocol TypeId, the variable starting from the m is kind ns3 convention to use it
	Ptr<Socket>     m_listningSocket;       // Listening socket
	std::list<Ptr<Socket> > m_acceptedsocketList; //the accepted sockets
	Ptr<Socket> sendingSocket;		// Socket for sending data
	Ptr<Socket> sharedSocket;
	uint32_t        m_totalbytesRecv;      // Total bytes received

	Ptr<RandomWaypointMobilityModel> m_mobility;	// mobility model of this node


	Ptr< WifiMacQueue > macQueue;		// pointer to the MAC queue

	std::queue<samples_t> samples_queue;
	std::queue<svm_node *> testfeatures_queue;
	/////////////////////////////////////////////////LWPR required variables initialization//////////////////////////////



//	LWPR_Object lwpr_model(23,1);

	double lwpr_mse_T;
	double lwpr_mse_P;
	int lwpr_counter;

	double lwpr_x[23];
	double lwpr_y;

	LWPR_Model lwpr_model;
//	doubleVec lwpr_x(23);
//	doubleVec lwpr_y(1);


	//actually am puting the samples in structure array, because otherwize its not working
	samples_t samples;
	samples_t ssamples; //shared samples array


	///////////////////////these four parameteres we need in our helper ///////////////////////////////
	Ipv4InterfaceContainer interfaces;	// known interfaces
	NodeContainer m_nodeContainer;	// known interfaces
	Ptr<Node> m_actualNodes;
	//	Ipv4Address m_source;

	//static std::string modelName;
	std::string outFilePath;		// path to the above file (for writing samples)
	bool mutex;
	// model parameters
	struct svm_model* model;
	struct svm_node *testfeatures;

	//struct svm_node *sharedfeatures;
	//	  struct svm_problem *problems;
	//	  struct svm_parameter *parameters;
	int nooffeatures;
	int svm_type;
	int nr_class;

	uint8_t check;
	double errSum;
	uint32_t errCounter;
	//////////////////////////////////////////////////////////////////////////////////////////////////
	//	static uint32_t m_sentbytes1;
	//	static uint32_t m_sentbytes2;
	//	static uint32_t m_sentbytes;

	std::map<Address ,uint32_t> m_Mapp;
	/* Bit-array for node status (neighborhod) */
	uint32_t my_status[5];
	/* Summary message counter */
	uint32_t counter;

	//shared samples with other node for colaborative learning
//	shared_samples_t *shared_samples;
//	shared_samples_t *recv_sharedsamples;
	double writefeatures[24];

	int ret;

	double comparing[24];

	double oldSamples[24];
	double recentSamples[24];

	online_rate_entry_t m_estimated_rate;
	/////End of Eduardo definitions///.....

	online_neighbor_entry_t Neighbors[MAX_NODE_ID];
	/* Sampling data table */
	online_sample_t Samples[MAX_NODE_ID];
	/* Last time node state was updated */
	uint64_t last_update;

	uint8_t oldRate;
	uint32_t oldEstRate;
	uint32_t oldBytesSent;
	//	UniformVariable uniformvariable(0.0,0.4);

	uint8_t hasEMA;
	uint8_t hasEWMA;
	float em_interval_sum;
	double interval;
	double starttime;

	summary_msg_t *summary;
	/* Sent packets counter since last summary */
	uint32_t sent_traffic_counter;
	/* Number of samples collected by the node*/
	uint32_t n_samples;
	uint32_t ss_counter; //counter for shared samples, how many we received so far

	uint32_t myid;

	EventId event_estimator;
	EventId m_sendEvent;    		// Eventid of pending "send packet" event
	EventId summaryTimer;
	EventId m_statusTimer;
	EventId ssampleEvent;

	Ipv4Address destination;
	Ipv4Address sharedDestination;

	double sent_time_counter;
	uint32_t m_summary_counter;
	//	static uint32_t m_my_status;
	uint64_t m_last_update;

	uint32_t m_totalNode;
	uint32_t m_individual_counter;

	double last_sent_time;
	//Time interval;
	double m_Simulationtime;

	FILE* fp;				// file for writing samples

	char Nodeid[10];
	std::string saved_sample;
	std::string model_name;
	std::string read_model;
	std::string predict_error;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
} //END OF namespace ns3

#endif /* MAPPLE_ONLINE_AAP_H_ */
