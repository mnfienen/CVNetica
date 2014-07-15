'''
Netica Constants taken from Netica.h --> Version 5.04
These should be updated with new versions of Netica.h to ensure compatability

Each Netica function for which specific constants are defined has been given
a class using the function name and appended with _const. 
General constants for which no function is identified are in the class
netica_const.
'''


class netica_const:
    REAL_VALUE = -25
    STATE_VALUE = -20
    GAUSSIAN_VALUE = -15
    INTERVAL_VALUE = -10
    STATE_NOT_VALUE = -7
    LIKELIHOOD_VALUE=-6
    NO_VALUE = -3
    EVERY_STATE = -5
    IMPOSS_STATE=-4
    UNDEF_STATE=-3
    FIRST_CASE = -15
    NEXT_CASE=-14
    NO_MORE_CASES=-13
    ENTROPY_SENSV = 0x02
    REAL_SENSV = 0x04
    VARIANCE_SENSV = 0x100
    VARIANCE_OF_REAL_SENSV = 0x104
    NEGATIVE_FINDING = -7
    LIKELIHOOD_FINDING=-6
    NO_FINDING = -3
    NO_VISUAL_INFO=0
    NO_WINDOW=0x10
    MINIMIZED_WINDOW=0x30
    REGULAR_WINDOW=0x70
    BELIEF_UPDATE = 0x100
    HALT_CALLBACK_RESULT = -1003
    ALL_THREADS = 0x20
    LAST_ENTRY = -10
    QUERY_ns = -1

class checking_nc_const:
    NO_CHECK=1
    QUICK_CHECK=2
    REGULAR_CHECK=3
    COMPLETE_CHEC=4
    QUERY_CHECK=-1

class errseverity_ns_const:
    NOTHING_ERR=1
    REPORT_ERR=2
    NOTICE_ERR=3
    WARNING_ERR=4
    ERROR_ERR=5
    XXX_ERR=6
    
class errcond_ns_const:
    OUT_OF_MEMORY_CND=0x08
    USER_ABORTED_CND=0x20
    FROM_WRAPPER_CND=0x40
    FROM_DEVELOPER_CND=0x80
    INCONS_FINDING_CND=0x200

class eventtype_ns_const:
    CREATE_EVENT=0x01
    DUPLICATE_EVENT=0x02
    REMOVE_EVENT=0x04

class nodetype_bn_const:
    CONTINUOUS_TYPE=1
    DISCRETE_TYPE=2
    TEXT_TYPE=3

class nodkind_bn_const:
    NATURE_NODE=1
    CONSTANT_NODE=2
    DECISION_NODE=3
    UTILITY_NODE=4
    DISCONNECTED_NODE=5
    ADVERSARY_NODE=6


class sampling_bn_const:
    DEFAULT_SAMPLING, JOIN_TREE_SAMPLING, FORWARD_SAMPLING = range(3)
    
class learn_method_bn_const:
    COUNTING_LEARNING=1
    EM_LEARNING=3
    GRADIENT_DESCENT_LEARNING=4

