#-----------------------------------------------------------------
# plyparser.py
#
# PLYParser class and other utilites for simplifying programming
# parsers with PLY
#
# Copyright (C) 2008-2011, Eli Bendersky
# License: BSD
#-----------------------------------------------------------------
import os
import re
import sys
import ply.lex
import logging
from ply.lex import TOKEN
import ply.yacc

class Coord(object):
    """ Coordinates of a syntactic element. Consists of:
            - File name
            - Line number
            - (optional) column number, for the Lexer
    """
    def __init__(self, file, line, column=None):
        self.file = file
        self.line = line
        self.column = column

    def __str__(self):
        str = "%s:%s" % (self.file, self.line)
        if self.column: str += ":%s" % self.column
        return str

class ParseError(Exception):
    pass

class PLYParser(object):
    def _create_opt_rule(self, rulename):
        """ Given a rule name, creates an optional ply.yacc rule
            for it. The name of the optional rule is
            <rulename>_opt
        """
        optname = rulename + '_opt'

        def optrule(self, p):
            p[0] = p[1]

        optrule.__doc__ = '%s : empty\n| %s' % (optname, rulename)
        optrule.__name__ = 'p_%s' % optname
        setattr(self.__class__, optrule.__name__, optrule)

    def _coord(self, lineno, column=None):
        return Coord(
                file=self.clex.filename,
                line=lineno,
                column=column)

    def _parse_error(self, msg, coord):
        raise ParseError("%s: %s" % (coord, msg))

#-----------------------------------------------------------------
# pycparser: c_lexer.py
#
# CLexer class: lexer for the C language
#
# Copyright (C) 2008-2011, Eli Bendersky
# License: BSD
#-----------------------------------------------------------------

class AbaqusLexer(object):
    """ A lexer for the C language. After building it, set the
        input text with input(), and call token() to get new
        tokens.

        The public attribute filename can be set to an initial
        filaneme, but the lexer will update it upon #line
        directives.
    """
    def __init__(self, error_func):
        """ Create a new Lexer.

            error_func:
                An error function. Will be called with an error
                message, line and column as arguments, in case of
                an error during lexing.

        """
        self.error_func = error_func
        self.filename = ''

    def build(self, **kwargs):
        """ Builds the lexer from the specification. Must be
            called after the lexer object is created.

            This method exists separately, because the PLY
            manual warns against calling lex.lex inside
            __init__
        """
        self.lexer = ply.lex.lex(object=self, **kwargs)

    def reset_lineno(self):
        """ Resets the internal line number counter of the lexer.
        """
        self.lexer.lineno = 1

    def input(self, text):
        self.lexer.input(text)

    def token(self):
        g = self.lexer.token()
        return g

    ######################--   PRIVATE   --######################

    ##
    ## Internal auxiliary methods
    ##
    def _error(self, msg, token):
        location = self._make_tok_location(token)
        self.error_func(msg, location[0], location[1])
        self.lexer.skip(1)

    def _find_tok_column(self, token):
        i = token.lexpos
        while i > 0:
            if self.lexer.lexdata[i] == '\n': break
            i -= 1
        return (token.lexpos - i) + 1

    def _make_tok_location(self, token):
        return (token.lineno, self._find_tok_column(token))

    ##
    ## All the tokens recognized by the lexer
    ##
    tokens = (
        # Keyword
        'KEYWORD',

        # Identifiers
        'ID',
        'PARAM',

        # constants
        'INT_CONST_DEC',
        'FLOAT_CONST',
        'CHAR_CONST',
        'WCHAR_CONST',

        # String literals
        'STRING_LITERAL',
        'WSTRING_LITERAL',

        # Assignment
        'EQUALS',

        # Delimeters
        'COMMA', 'PERIOD',          # . ,

        'LASTTOKENONLINE',
    )

    ##
    ## Regexes for use in tokens
    ##
    ##

    # valid C identifiers (K&R2: A.2.3)
    identifier = r'[a-zA-Z_][0-9a-zA-Z_. ]*'
    abaqus_keyword = r'\*[a-zA-Z][0-9a-zA-Z_ \t]*'
    abaqus_identifier = r'[a-zA-Z_][0-9a-zA-Z_. \-\:\/]*'
    identifier = r'[a-zA-Z_][0-9a-zA-Z_. \-\:\/]*'

    # integer constants (K&R2: A.2.5.1)
    integer_suffix_opt = r'(u?ll|U?LL|([uU][lL])|([lL][uU])|[uU]|[lL])?'
    decimal_constant = '(0'+integer_suffix_opt+')|([1-9][0-9]*'+integer_suffix_opt+')'
    octal_constant = '0[0-7]*'+integer_suffix_opt
    hex_constant = '0[xX][0-9a-fA-F]+'+integer_suffix_opt

    bad_octal_constant = '0[0-7]*[89]'

    # character constants (K&R2: A.2.5.2)
    # Note: a-zA-Z and '.-~^_!=&;,' are allowed as escape chars to support #line
    # directives with Windows paths as filenames (..\..\dir\file)
    #
    simple_escape = r"""([a-zA-Z._~!=&\^\-\\?'"])"""
    octal_escape = r"""([0-7]{1,3})"""
    hex_escape = r"""(x[0-9a-fA-F]+)"""
    bad_escape = r"""([\\][^a-zA-Z._~^!=&\^\-\\?'"x0-7])"""

    escape_sequence = r"""(\\("""+simple_escape+'|'+octal_escape+'|'+hex_escape+'))'
    cconst_char = r"""([^'\\\n]|"""+escape_sequence+')'
    char_const = "'"+cconst_char+"'"
    wchar_const = 'L'+char_const
    unmatched_quote = "('"+cconst_char+"*\\n)|('"+cconst_char+"*$)"
    bad_char_const = r"""('"""+cconst_char+"""[^'\n]+')|('')|('"""+bad_escape+r"""[^'\n]*')"""

    # string literals (K&R2: A.2.6)
    string_char = r"""([^"\\\n]|"""+escape_sequence+')'
    string_literal = '"'+string_char+'*"'
    wstring_literal = 'L'+string_literal
    bad_string_literal = '"'+string_char+'*'+bad_escape+string_char+'*"'

    # floating constants (K&R2: A.2.5.3)
    exponent_part = r"""([eE][-+]?[0-9]+)"""
    fractional_constant = r"""([\+\-]*[0-9]*\.[0-9]+)|([\+\-]*[0-9]+\.)"""
    floating_constant = '(((('+fractional_constant+')'+exponent_part+'?)|([0-9]+'+exponent_part+'))[FfLl]?)'

    ##
    ## Lexer states
    ##
    states = (
              #
              ('keywordstate', 'inclusive'),
              ('datalinestate', 'inclusive'),
             )

    @TOKEN(abaqus_keyword)
    def t_KEYWORD(self, t):
        t.lexer.push_state('keywordstate')
        t.value = t.value[1:]
        return t

    @TOKEN(abaqus_identifier)
    def t_keywordstate_PARAM(self, t):
        return t

    def t_keywordstate_CONTINUE(self, t):
        r',\n'
        t.value = ','
        t.type = 'COMMA'
        return t

    def t_keywordstate_NEWLINE(self, t):
        r'\n'
        t.lexer.pop_state()
        t.lexer.push_state('datalinestate')
        t.lexer.lineno += t.value.count("\n")

    t_keywordstate_ignore = ' \t'

    def t_keywordstate_error(self, t):
        msg = 'INVALID KEYWORD'
        self._error(msg, t)

    def t_datalinestate_LASTTOKENONLINE(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")
        return t

    ##
    ## Rules for the normal state
    ##
    t_ignore = ' \t'


    def t_COMMENT(self, t):
        r'[ \t]*\*\*.*\n'
        t.lexer.lineno += t.value.count("\n")

    # Newlines

    def t_NEWLINEALONE(self, t):
        r'\n+'
        t.lexer.lineno += t.value.count("\n")

    # Assignment operators
    t_EQUALS            = r'='

    # Delimeters
    t_COMMA             = r','
    t_PERIOD            = r'\.'

    t_STRING_LITERAL    = string_literal

    # The following floating and integer constants are defined as
    # functions to impose a strict order (otherwise, decimal
    # is placed before the others because its regex is longer,
    # and this is bad)
    #
    @TOKEN(floating_constant)
    def t_FLOAT_CONST(self, t):
        return t

    #@TOKEN(hex_constant)
    #def t_INT_CONST_HEX(self, t):
    #    return t

    #@TOKEN(bad_octal_constant)
    #def t_BAD_CONST_OCT(self, t):
    #    msg = "Invalid octal constant"
    #    self._error(msg, t)

    #@TOKEN(octal_constant)
    #def t_INT_CONST_OCT(self, t):
    #    return t

    @TOKEN(decimal_constant)
    def t_INT_CONST_DEC(self, t):
        return t

    # Must come before bad_char_const, to prevent it from
    # catching valid char constants as invalid
    #
    @TOKEN(char_const)
    def t_CHAR_CONST(self, t):
        return t

    @TOKEN(wchar_const)
    def t_WCHAR_CONST(self, t):
        return t

    @TOKEN(unmatched_quote)
    def t_UNMATCHED_QUOTE(self, t):
        msg = "UNMATCHED '"
        self._error(msg, t)

    @TOKEN(bad_char_const)
    def t_BAD_CHAR_CONST(self, t):
        msg = "INVALID CHAR CONSTANT %s" % t.value
        self._error(msg, t)

    @TOKEN(wstring_literal)
    def t_WSTRING_LITERAL(self, t):
        return t

    # unmatched string literals are caught by the preprocessor

    @TOKEN(bad_string_literal)
    def t_BAD_STRING_LITERAL(self, t):
        msg = "STRING CONTAINS INVALID ESCAPE CODE"
        self._error(msg, t)

    @TOKEN(identifier)
    def t_ID(self, t):
        t.type = "ID"
        return t

    def t_error(self, t):
        msg = 'ILLEGAL CHARACTER %s' % repr(t.value[0])
        self._error(msg, t)

def list_append(lst, item):
    lst.append(item)
    return lst

def chunker(seq, size=8):
    return (seq[pos:pos + size] for pos in range(0, len(seq), size))

def cat(key, sep):
    return sep.join([s.strip() for s in key.split(' ') if s.split()])

class Parameter(object):
    def __init__(self, name, value = None):
        self.name = name.strip()
        self.value = value
        self.key = cat(self.name, '_').lower()

class Parameters(object):
    def __init__(self, params):
        self.data = params or []
    def __iter__(self):
        return iter(self.data)
    def __getitem__(self, i):
        if isinstance(i, int):
            return self.data[i]
        return self.get(i)
    def __contains__(self, key):
        return key.lower() in [p.key for p in self.data]
    def __bool__(self):
        return bool(self.data)
    def get(self, key):
        kw = key.lower()
        for p in self.data:
            if kw == p.key:
                return p.value
    def todict(self):
        return dict([(p.key, p.value) for p in self.data])
    def toinp(self):
        return ', '.join(['{0}={1}'.format(p.key, p.value) for p in self.data])
    def index(self, key):
        kw = key.lower()
        for (i, p) in enumerate(self.data):
            if kw == p.key:
                return i
        raise IndexError(key)

class Keyword(object):
    def __init__(self, keyword, params=None, data=None):
        self.keyword = keyword.strip()
        self.params = Parameters(params)
        self.data = data
        self.key = cat(self.keyword, '_').lower()

    def __str__(self):
        kwd_str_list = ['Keyword:{0}'.format(self.keyword),]
        kwd_str_list.append('\n')
        for param in self.params:
            kwd_str_list.append('Parameter:{0}={1}'.format(param.name, param.value))
            kwd_str_list.append('\n')
        kwd_str_list.append('Data:\n')
        if self.data is not None:
            for data in self.data:
                kwd_str_list.append('{0},'.format(data))
                kwd_str_list.append('\n')
        return ''.join([item for item in kwd_str_list])

    def __repr__(self):
        return str(self)

    def toinp(self):
        o = ['*{0}'.format(self.keyword.title())]
        if self.params:
            o[0] += ', ' + self.params.toinp()
        if self.data is not None:
            for data in self.data:
                chunks = chunker(data)
                for chunk in chunker(data):
                    o.append(' ' + ', '.join(chunk) + ',')
                o[-1] = o[-1][:-1]
        return '\n'.join(o)

keywords = []

class AbaqusParser(PLYParser):
    def __init__(
            self,
            lex_optimize=True,
            lextab='pycparser.lextab',
            yacc_optimize=True,
            yacctab='pycparser.yacctab',
            yacc_debug=False):
        """ Create a new AbaqusParser.

            Some arguments for controlling the debug/optimization
            level of the parser are provided. The defaults are
            tuned for release/performance mode.
            The simple rules for using them are:
            *) When tweaking AbaqusParser/CLexer, set these to False
            *) When releasing a stable parser, set to True

            lex_optimize:
                Set to False when you're modifying the lexer.
                Otherwise, changes in the lexer won't be used, if
                some lextab.py file exists.
                When releasing with a stable lexer, set to True
                to save the re-generation of the lexer table on
                each run.

            lextab:
                Points to the lex table that's used for optimized
                mode. Only if you're modifying the lexer and want
                some tests to avoid re-generating the table, make
                this point to a local lex table file (that's been
                earlier generated with lex_optimize=True)

            yacc_optimize:
                Set to False when you're modifying the parser.
                Otherwise, changes in the parser won't be used, if
                some parsetab.py file exists.
                When releasing with a stable parser, set to True
                to save the re-generation of the parser table on
                each run.

            yacctab:
                Points to the yacc table that's used for optimized
                mode. Only if you're modifying the parser, make
                this point to a local yacc table file

            yacc_debug:
                Generate a parser.out file that explains how yacc
                built the parsing table from the grammar.
        """
        self.clex = AbaqusLexer(error_func=self._lex_error_func)

        self.clex.build(
            optimize=lex_optimize,
            lextab=lextab)
        self.tokens = self.clex.tokens

        self.cparser = ply.yacc.yacc(
            module=self,
            start='keyword_list',
            debug=yacc_debug,
            optimize=yacc_optimize,
            tabmodule=yacctab)


    def parse(self, text, filename='', debuglevel=0):
        """ Parses Abaqus input files and returns an AST.

            text:
                A string containing the Abaqus input file

            filename:
                Name of the file being parsed (for meaningful
                error messages)

            debuglevel:
                Debug level to yacc
        """
        self.clex.filename = filename
        self.clex.reset_lineno()
        try:
            text = text.decode()
        except AttributeError:
            pass
        t = self.cparser.parse(text, lexer=self.clex, debug=debuglevel)
        return t

    ######################--   PRIVATE   --######################


    def _lex_error_func(self, msg, line, column):
        self._parse_error(msg, self._coord(line, column))

    ##
    ## Grammar productions
    ##

    def p_keyword_list(self, p):
        '''
        keyword_list : keyword_list keyword
        '''
        p[0] = p[1] + [p[2]]

    def p_keyword(self, p):
        '''
        keyword_list : keyword
        '''
        p[0] = [p[1]]

    def p_single_keyword(self, p):
        '''
        keyword : KEYWORD
                | KEYWORD data_lines
                | KEYWORD COMMA param_list
                | KEYWORD COMMA param_list data_lines
        '''
        if len(p) == 2:
            # KEYWORD
            p[0] = Keyword(p[1])
        elif len(p) == 3:
            # KEYWORD data_list
            p[0] = Keyword(p[1], data = p[2])
        elif len(p) == 4:
            # KEYWORD COMMA param_list
            p[0] = Keyword(p[1], params = p[3])
        elif len(p) == 5:
            # KEYWORD COMMA param_list data_list
            p[0] = Keyword(p[1], params = p[3], data = p[4])
        else:
            # Error?
            pass

    def p_param_list(self, p):
        '''param_list : param_list COMMA param'''
        p[0] = p[1] + [p[3]]

    def p_param(self, p):
        '''param_list : param'''
        p[0] = [p[1]]

    def p_single_param(self, p):
        '''
        param : PARAM
              | PARAM EQUALS PARAM
              | PARAM EQUALS FLOAT_CONST
              | PARAM EQUALS INT_CONST_DEC
        '''
        if len(p) == 2:
            p[0] = Parameter(p[1])
        elif len(p) == 4:
            p[0] = Parameter(p[1], value = p[3])

    def p_data_lines_list(self, p):
        '''
        data_lines : data_lines data_line
        '''
        p[0] = list_append(p[1],p[2])

    def p_data_lines(self, p):
        '''
        data_lines : data_line
        '''
        p[0] = [p[1]]

    def p_data_line(self, p):
        '''
        data_line : data_list LASTTOKENONLINE
                  | data_list COMMA LASTTOKENONLINE
        '''
        p[0] = p[1]

    def p_data_list(self, p):
        '''
        data_list : data_list COMMA data
                  | data_list COMMA LASTTOKENONLINE data
                  | data_list data
        '''
        if len(p) == 3:
            p[0] = p[1] + [p[2]]
        elif len(p) == 4:
            p[0] = p[1] + [p[3]]
        elif len(p) == 5:
            p[0] = p[1] + [p[4]]

    def p_data(self, p):
        '''
        data_list : data
        '''
        p[0] = [p[1]]

    def p_single_data(self, p):
        '''
        data : ID
             | INT_CONST_DEC
             | FLOAT_CONST
        '''
        p[0] = p[1]

    def p_error(self, p):
        if p:
            self._parse_error(
                'BEFORE: %s' % p.value,
                self._coord(p.lineno))
        else:
            self._parse_error('AT END OF INPUT', '')

if __name__ == "__main__":
    text = '''
    *heading
    word1 word2
    line2
    ** comment
    *keyword,singleparam
    *KEYword,
    param=continue
    *KEYword,param=coffee,param=1.0,param=3,param=4.0e-3
    *node,nset=all_nodes
    1,1.0,1.0e-5,1.0E+6
    2,1.0,1.0e-5,1.0E+6
    3,1.0,1.0e-5,1.0E+6
    4,1.0,1.0e-5,1.0E+6
    *element,type=c3d4,elset=foo:1
    1,1,2,3,4
    2,3,4,5,6
    *element,type=c3d15,elset=foo:2
    1,1,2,3,4,5,6,7,
    8,9,10,11,12,13,14,15
    *end
    '''
    def errfoo(msg, a, b):
        print(msg + "\n")
        sys.exit()

    clex = AbaqusLexer(errfoo)
    clex.build(reflags=re.IGNORECASE)
    clex.input(text)

    while 1:
        tok = clex.token()
        if not tok: break

        #~ print type(tok)
        print([tok.value, tok.type, tok.lineno, clex.filename, tok.lexpos])
