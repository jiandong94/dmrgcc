#include "util/input.h"

void InputFile::open()
{
    close();
    file_.open(filename_.c_str());
    if(!file_)
        error("can't open input file");
    opened_ = true;
}

void InputFile::close()
{
    if(opened_)
    {
        file_.clear();
        file_.close();
    }
    opened_ = false;
}

//----------------------------------------------------------------------------------------

void eatwhite(istream& is)
{
    char c;
    while(is.get(c))
    {
        if(!isspace(c))
        {
            is.putback(c);
            break;
        }
    }
}

void matchbracket(istream& is)
{
    int nest = 1;
    char c;
    while(is >> c)
    {
        if(c == '{')
            nest++;
        else if(c == '}')
            nest--;
        if(nest == 0) break;
    }
    if(nest != 0)
        error("unterminated bracket");
}

int gettoken(istream& is, string& s)
{
    char c[512], ch;
    eatwhite(is);
    for(int i=0;i<512;++i)
    {
        if(!is.get(ch))
            return 0;
        if(isalpha(ch) || ch == '_' || (i!=0 && isdigit(ch)))
            c[i] = ch;
        else
        {
            if(i == 0)
            {
                c[0] = ch;
                c[1] = '\0';
                s = c;
                return 1;
            }
            is.putback(ch);
            c[i] = '\0';
            break;
        }
    }
    s = c;
    return 1;
}

//----------------------------------------------------------------------------------------

InputGroup::InputGroup(string file_name, string group_name)
: parent_(0), name_(group_name)
{
    input_file_ = new InputFile(file_name);
}

InputGroup::InputGroup(InputGroup& par, string name)
: parent_(&par), name_(name)
{
    input_file_ = &(*par.input_file_);
}

InputGroup::~InputGroup() { }

int InputGroup::GotoGroup()
{
    if(parent_ != 0)
        parent_->GotoGroup();
    else 
        input_file_->open();
    eatwhite(input_file_->file());
    while(1)
    {
        string s;
        if(!gettoken(input_file_->file(), s))
            return 0;
        if(s == name_)
        {
            eatwhite(input_file_->file());
            char c;
            input_file_->file().get(c);
            if(c != '{')
                error("bracket does not follow group name");
            eatwhite(input_file_->file());
            return 1;
        }
        else
        {
            if(s[0] == '{')
            {
                matchbracket(input_file_->file());
            }
            while(1)
            {
                char c;
                input_file_->file().get(c);
                if(c == '{')
                    matchbracket(input_file_->file());
                if(c == '\n')
                {
                    eatwhite(input_file_->file());
                    break;
                }
            }

        }
    }
}


int InputGroup::GotoToken(string s)
{
    if(!GotoGroup())
        return 0;
    while(1)
    {
        string t;
        if(!gettoken(input_file_->file(), t))
        {
            return 0;
        }
        if(t[0] == '}') return 0;
        if(t == s)
        {
            eatwhite(input_file_->file());
            char c;
            if(!(input_file_->file().get(c))) return 0;
            if(c != '=') return 0;
            eatwhite(input_file_->file());
            return 1;
        }
        if(t[0] == '{')
        {
            matchbracket(input_file_->file());
        }
        while(1)
        {
            char c;
            input_file_->file().get(c);
            if(c == '{')
                matchbracket(input_file_->file());
            if(c == '\n')
            {
                eatwhite(input_file_->file());
                break;
            }

        }
    }
}

void InputGroup::SkipLine()
{
    char c = '\0';
    while(c != '\n')
        input_file_->file().get(c);
    eatwhite(input_file_->file());
}

int InputGroup::GetInt(string s, int& res, bool hasdf, int df)
{
    if(!GotoToken(s) || !(input_file_->file() >> res))
    {
        if(hasdf) printf("Def %s.%s = %d\n", name_.c_str(), s.c_str(), df);
        return 0;
    }
    printf("Got %s.%s = %d\n", name_.c_str(), s.c_str(), res);
    return 1;
}

int InputGroup::GetDouble(string s, double& res, bool hasdf, double df)
{
    if(!GotoToken(s) || !(input_file_->file() >> res))
    {
        if(hasdf) printf("Def %s.%s = %.2f\n", name_.c_str(), s.c_str(), df);
        return 0;
    }
    printf("Got %s.%s = %.2f\n", name_.c_str(), s.c_str(), res);
    return 1;
}

int InputGroup::GetString(string s, string& res, bool hasdf, string df)
{
    if(!GotoToken(s) || !(input_file_->file() >> res))
    {
        if(hasdf) printf("Def %s.%s = %s\n", name_.c_str(), s.c_str(), df.c_str());
        return 0;
    }
    printf("Got %s.%s = %s\n", name_.c_str(), s.c_str(), res.c_str());
    return 1;
}

int InputGroup::GetBool(string s, bool& res, bool hasdf, bool df)
{
    string tmp_res;
    if(!GotoToken(s) || !(input_file_->file() >> tmp_res))
    {
        if(hasdf) printf("Def %s.%s = %d\n", name_.c_str(), s.c_str(), df);
        return 0;
    }
    if(tmp_res == "true" || tmp_res == "1")
        res = 1;
    else if(tmp_res == "false" || tmp_res == "0")
        res = 0;
    printf("Got %s.%s = %d\n", name_.c_str(), s.c_str(), res);
    return 1;
}

int InputGroup::getInt(string s, int def)
{
    int res = 0;
    int got = GetInt(s, res, true, def);
    if(!got) return def;
    else return res;
}

double InputGroup::getDouble(string s, double def)
{
    double res = 0.0;
    int got = GetDouble(s, res, true, def);
    if(!got) return def;
    else return res;
}

string InputGroup::getString(string s, string def)
{
    string res;
    int got = GetString(s, res, true, def);
    if(!got) return def;
    else return res;
}

bool InputGroup::getBool(string s, bool def)
{
    bool res = false;
    int got = GetBool(s, res, true, def);
    if(!got) return def;
    else return res;
}

int InputGroup::getInt(string s)
{
    int res = 0;
    if(!GetInt(s, res))
        error("there is no item: " + s);
    return res;
}

double InputGroup::getDouble(string s)
{
    double res = 0.0;
    if(!GetDouble(s, res))
        error("there is no item: " + s);
    return res;
}

string InputGroup::getString(string s)
{
    string res;
    if(!GetString(s, res))
        error("there is no item: " + s);
    return res;
}

bool InputGroup::getBool(string s)
{
    bool res = false;
    if(!GetBool(s, res))
        error("there is no item: " + s);
    return res;
}
