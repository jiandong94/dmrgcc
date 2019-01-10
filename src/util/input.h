#ifndef DMRGCC_UTIL_INPUT_H_
#define DMRGCC_UTIL_INPUT_H_

#include "util/general.h"
class InputFile
{
    private:

    string filename_;

    ifstream file_;

    bool opened_;

    public:

    InputFile() : opened_(false) { }

    InputFile(string fname) : filename_(fname), opened_(false) { }

    void open();

    void close();

    string const& filename() 
           const { return filename_; }

    ifstream& file()
              { return file_; }

    ifstream const& file()
             const { return file_; }

    bool opened()
         const { return opened_; }
};

class InputGroup
{
    private:

    InputFile* input_file_;

    InputGroup* parent_;

    string name_;

    public:

    InputGroup() {}

    InputGroup(string file_name, string group_name);

    InputGroup(InputGroup& parent, string name);

    ~InputGroup();

    const string& name() const { return name_; } 
    
    ifstream& file() const {return input_file_->file();}
   
    //
    int GotoGroup();

    int GotoToken(string s);

    void SkipLine();

    // return value
    int GetInt(string s, int& i, bool hasdf = false, int df = 0);
    int GetDouble(string s, double& d, bool hasdf = false, double df = 0.0);
    int GetString(string s, string& t, bool hasdf = false, string df = "");
    int GetBool(string s, bool& b, bool hasdf = false, bool df = false);

    int getInt(string s, int default_);
    double getDouble(string s, double default_);
    string getString(string s, string default_);
    bool getBool(string s, bool default_);

    int getInt(string s);
    double getDouble(string s);
    string getString(string s);
    bool getBool(string s);



};

#endif // DMRGCC_UTIL_INPUT_H_
