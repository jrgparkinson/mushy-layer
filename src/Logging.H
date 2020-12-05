/** Logging macros
 * Different levels of logging depending on s_verbosity (main.verbosity):
 * ERROR: s_verbosity >=-999
 * WARNING: s_verbosity >= 0
 * INFO: s_verbosity >= 1
 * DEBUG: s_verbosity >= 2
 * FUNCTION: s_verbosity >= 3 (logs when we enter various functions)
 * DEBUG EXTRA: s_verbosity >= 4 even more debug logging, lots and lots of output
**/
#ifndef LOG

#define __FILENAME__ (strrchr(__FILE__, '/') ? strrchr(__FILE__, '/') + 1 : __FILE__)
#define LOG_BASE(msg, prefix, verb_threshold, verb, end) do{  \
            if (verb >= verb_threshold) { \
                pout() << "[" << __FILENAME__ << ":" << __LINE__ << "] " << prefix << msg; \
                if (end) { pout() << endl; } \
            } \
            } while( false )

// For use in classes without s_verbosity defined
#define LOG(msg) do{  LOG_BASE(msg, "", 0, 1, true); } while( false )
#define LOG_NOEND(msg) do{  LOG_BASE(msg, "", 0, 1, false); } while( false )

// For logging in AMRLevelMushyLayer
#define LOG_MAIN_CLASS(msg, prefix, verb) do{  LOG_BASE(msg, prefix, verb, s_verbosity,true);  } while( false ) 
#define LOG_DEBUG_EXTRA(msg) do{  LOG_MAIN_CLASS(msg, "", 4);  } while( false )
#define LOG_DEBUG(msg) do{  LOG_MAIN_CLASS(msg, "", 2);  } while( false )
#define LOG_INFO(msg) do{ LOG_BASE(msg, "", 1, s_verbosity, true);  } while( false )
#define LOG_WARNING(msg) do{ LOG_BASE(msg, "WARNING: ", 0, s_verbosity, true); } while( false )
#define LOG_ERROR(msg) do{ LOG_BASE(msg, "ERROR: ", -999, s_verbosity, true); } while( false )

#define LOG_FUNCTION(function_desc, entry_exit) do{ LOG_BASE(function_desc << " (level " << m_level << ")", entry_exit, 3, s_verbosity, true); } while( false )
#define LOG_FUNCTION_ENTRY() do{ LOG_FUNCTION(__FUNCTION__, "Entering: "); } while( false )
#define LOG_FUNCTION_EXIT() do{ LOG_FUNCTION(__FUNCTION__, "Leaving: "); } while( false )

#endif