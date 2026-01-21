#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <map>
#include <span>
#include <variant>
#include <algorithm>
#include <chrono>

using Index = std::vector<int>;

struct Token
{
	uint32_t word;
	uint32_t c5;
	uint32_t lemma;
	uint32_t pos;
};

struct Literal
{
	std::string attribute;
	uint32_t value; // CHANGED
	bool is_equality;
};

struct Corpus
{
	std::vector<Token> tokens;
	std::vector<int> sentences;
	std::vector<std::string> index2string;
	std::map<std::string, uint32_t> string2index;
	Index word_index;  // NEW
	Index c5_index;	   // NEW
	Index lemma_index; // NEW
	Index pos_index;   // NEW
};

struct IndexSet
{
	std::span<const int> elems;
	int shift; // NEW

	IndexSet()
	    : elems(std::span<const int>()), 
	      shift(0)
	{
	}

	IndexSet(const std::span<const int> &elems, int shift)
	    : elems(elems), 
	      shift(shift)
	{
	} 
};

struct DenseSet
{
	int first;
	int last;

	DenseSet(int start, int end)
	    : first(start), last(end)
	{
	}
};

struct ExplicitSet
{
	std::vector<int> elems; 

	ExplicitSet() : elems() {} 
};

struct MatchSet
{
	std::variant<DenseSet, IndexSet, ExplicitSet> set;
	bool complement;

	MatchSet()
	    : set(IndexSet()), 
	      complement(false)
	{
	} 
};

using Sentence = std::vector<Token>;
using Clause = std::vector<Literal>;
using Query = std::vector<Clause>;

typedef struct
{
	int sentence; // the matching sentence
	int pos;      // the first token in the match (indexing into sentence)
	int len;      // the number of consecutive tokens in the match (= # clauses in the query)
} Match;

std::vector<Match> match(const Corpus &corpus, const std::string &query_string);
Corpus load_corpus(const std::string &filename);
Query parse_query(const std::string &text, Corpus &corpus);
std::vector<Match> match(const Corpus &corpus, const Query &query);
void find_Highlights(bool &pos_h, bool &lemma_h, bool &word_h, bool &c5_h, const Clause clause, Token token);
void print_token(const Token token, const Clause clause, bool find_highlight, const Corpus &corpus);
void DisplayMatch(Corpus corpus, Match match, Query query);
Index build_index(const std::vector<Token> &tokens, uint32_t Token::*attribute);
void build_indices(Corpus &corpus);
IndexSet index_lookup(const Corpus &corpus, const std::string &attribute, uint32_t value);
std::vector<Match> match_single(const Corpus &corpus, const std::string &attr, const std::string &value);
MatchSet intersection(const MatchSet &A, const MatchSet &B);
ExplicitSet intersect(const ExplicitSet &A, const DenseSet &B);
ExplicitSet intersect(const ExplicitSet &A, const IndexSet &B);
ExplicitSet intersect(const DenseSet &A, const IndexSet &B);
ExplicitSet intersect(const ExplicitSet &A, const ExplicitSet &B);
ExplicitSet intersect(const IndexSet &A, const ExplicitSet &B);
ExplicitSet intersect(const DenseSet &A, const ExplicitSet &B);
ExplicitSet intersect(const IndexSet &A, const IndexSet &B);
DenseSet intersect(const DenseSet &A, const DenseSet &B);
ExplicitSet intersect(const IndexSet &A, const DenseSet &B);
ExplicitSet difference(const ExplicitSet &A, const ExplicitSet &B);
ExplicitSet difference(const IndexSet &A, const IndexSet &B);
ExplicitSet difference(const IndexSet &A, const ExplicitSet &B);
ExplicitSet difference(const ExplicitSet &A, const IndexSet &B);
DenseSet difference(const DenseSet &A, const DenseSet &B);
ExplicitSet difference(const DenseSet &A, const ExplicitSet &B);
ExplicitSet difference(const DenseSet &A, const IndexSet &B);
ExplicitSet difference(const ExplicitSet &A, const DenseSet &B);
ExplicitSet difference(const IndexSet &A, const DenseSet &B);
MatchSet match_set(const Corpus &corpus, const Literal &literal, int shift);
MatchSet match_set(const Corpus &corpus, const Clause &clause, int shift);
MatchSet match_set(const Corpus &corpus, const Query &query);
void add_match_if_valid(std::vector<Match> &matches, const Corpus &corpus, int set_index, int query_size);
std::vector<Match> match2(const Corpus &corpus, const Query &query);
void sort_matches(std::vector<MatchSet> &matches);
void split_matches(std::vector<MatchSet> &match_sets, std::vector<MatchSet> &d_sets, std::vector<MatchSet> &ie_sets);
MatchSet collapse_densets(std::vector<MatchSet> &d_sets);
ExplicitSet intersect_binary(const IndexSet &A, const IndexSet &B);
ExplicitSet intersect_binary(const IndexSet &A, const ExplicitSet &B);
ExplicitSet intersect_binary(const ExplicitSet &A, const ExplicitSet &B);
ExplicitSet difference_binary(const ExplicitSet &A, const ExplicitSet &B);
ExplicitSet difference_binary(const IndexSet &A, const IndexSet &B);
ExplicitSet difference_binary(const IndexSet &A, const ExplicitSet &B);
ExplicitSet difference_binary(const ExplicitSet &A, const IndexSet &B);

			

int main(void)
{
	std::cout<< "Loading corpus... Wait for instructions"<< std::endl;
	Corpus corpus = load_corpus("bnc-05M.csv");
	build_indices(corpus);
	Query query;

	do
	{
		std::cout << "Enter query, press enter to exit: ";
		std::string query_string;
		getline(std::cin, query_string);

		if (query_string.empty()) // if nothing entered, exit
		{
			std::cout << "No more queries. Exiting.... " << std::endl;
			break;
		}

		try
		{
			query = parse_query(query_string, corpus);
		}
		catch (const std::runtime_error &e)
		{
			std::cout << "Error: " << e.what() << std::endl;
			continue;
		}
	
		std::vector<Match> matches = match2(corpus, query);
		
		if (matches.empty())
		{
			std::cout << "No matches found." << std::endl;
		}
		else
		{
			for (size_t matches_index = 0; matches_index < matches.size() && matches_index < 10; matches_index++)
			{
				DisplayMatch(corpus, matches[matches_index], query);
			}
		}

	} while (1);
}

void print_token(const Token token, const Clause clause, bool find_highlight, const Corpus &corpus)
{
	bool lemma_h = false, pos_h = false, c5_h = false, word_h = false;

	if (find_highlight) // if we know we dont have a highlight on this token, skip finding highlights
		find_Highlights(pos_h, lemma_h, word_h, c5_h, clause, token);

	if (word_h)
		std::cout << "\033[1;31m" << corpus.index2string[token.word] << "\033[0m\t";
	else
		std::cout << corpus.index2string[token.word] << "\t";

	if (c5_h)
		std::cout << "\033[1;31m" << corpus.index2string[token.c5] << "\033[0m\t";
	else
		std::cout << corpus.index2string[token.c5] << "\t";

	if (lemma_h)
		std::cout << "\033[1;31m" << corpus.index2string[token.lemma] << "\033[0m\t";
	else
		std::cout << corpus.index2string[token.lemma] << "\t";

	if (pos_h)
		std::cout << "\033[1;31m" << corpus.index2string[token.pos] << "\033[0m\n";
	else
		std::cout << corpus.index2string[token.pos] << std::endl;
}

void find_Highlights(bool &pos_h, bool &lemma_h, bool &word_h, bool &c5_h, const Clause clause, Token token)
{
	for (size_t i = 0; i < clause.size(); i++)
	{
		if (clause[i].attribute == "pos")
		{
			if (clause[i].is_equality ? token.pos == clause[i].value : token.pos != clause[i].value)
				pos_h = true;
		}
		else if (clause[i].attribute == "c5")
		{
			if (clause[i].is_equality ? token.c5 == clause[i].value : token.c5 != clause[i].value)
				c5_h = true;
		}
		else if (clause[i].attribute == "lemma")
		{
			if (clause[i].is_equality ? token.lemma == clause[i].value : token.lemma != clause[i].value)
				lemma_h = true;
		}
		else if (clause[i].attribute == "word")
		{
			if (clause[i].is_equality ? token.word == clause[i].value : token.word != clause[i].value)
				word_h = true;
		}
	}
}

void DisplayMatch(Corpus corpus, Match match, Query query)
{
	int sentence_start = corpus.sentences[match.sentence];
	int sentence_end = corpus.sentences[match.sentence + 1];
	int clause_index = 0;

	int end_pos = match.pos + match.len + sentence_start;
	for (int sentence_index = sentence_start; sentence_index < sentence_end; sentence_index++)
	{
		Token current_token = corpus.tokens[sentence_index];
		bool highlight = (sentence_index >= (match.pos + sentence_start) && sentence_index < end_pos);

		if (highlight)
		{		
			print_token(current_token, query[clause_index], true, corpus);
			clause_index++;
		}
		else
			print_token(current_token, query[0], false, corpus); // send first clause since it always exists, wont be used either way.
	}
	std::cout << std::endl;
}
Corpus load_corpus(const std::string &filename)
{
	Corpus corpus;
	Token token;
	std::string line;
	std::vector<std::string> read_data(4);
	int token_index = 0; // Keeps track of token positions

	std::ifstream infile(filename);

	while (getline(infile, line))
	{
		if (line.empty())
		{
			if (token_index > 0)
			{
				// Mark the start of the next sentence
				corpus.sentences.push_back(token_index);
			}
			continue;
		}

		if (line[0] == '#') // Skip comment lines
			continue;

		std::stringstream ss(line);
		for (size_t i = 0; i < 4; i++)
			getline(ss, read_data[i], '\t');

		for (size_t j = 0; j < 4; j++)
		{
			if (corpus.string2index.find(read_data[j]) == corpus.string2index.end())
			{
				corpus.index2string.push_back(read_data[j]);
				corpus.string2index.insert(std::make_pair(read_data[j], corpus.index2string.size() - 1));
			}
			if (j == 0)
				token.word = corpus.string2index.find(read_data[j])->second;
			else if (j == 1)
				token.c5 = corpus.string2index.find(read_data[j])->second;

			else if (j == 2)
				token.lemma = corpus.string2index.find(read_data[j])->second;

			else if (j == 3)
				token.pos = corpus.string2index.find(read_data[j])->second;
		}

		corpus.tokens.push_back(token);
		token_index++;
	}

	if (!corpus.tokens.empty())
	{
		corpus.sentences.push_back(token_index);
	}

	return corpus;
}

Query parse_query(const std::string &text, Corpus &corpus)
{
	Literal literal;
	Query query;

	size_t stringPosition = 0;

	while (stringPosition < text.length())
	{
		Clause clause;

		while (stringPosition < text.length() && isspace(text[stringPosition]))
			stringPosition++;

		if (text[stringPosition] != '[')
			throw std::runtime_error("Syntax error: expected '[' at position " + std::to_string(stringPosition));

		stringPosition++;

		while (stringPosition < text.length() && text[stringPosition] != ']')
		{
			while (stringPosition < text.length() && isspace(text[stringPosition]))
				stringPosition++;

			int attributeStart = stringPosition;

			while (stringPosition < text.length() && isalnum(text[stringPosition]))
				stringPosition++;

			std::string tempAttribute = text.substr(attributeStart, stringPosition - attributeStart);

			if (tempAttribute != "word" && tempAttribute != "c5" && tempAttribute != "lemma" && tempAttribute != "pos")
				throw std::runtime_error("Syntax error: Attribute doesnt match any of the accepted values");

			literal.attribute = tempAttribute;
			if (stringPosition < text.length() && text[stringPosition] == '=')
			{
				literal.is_equality = true;
				stringPosition++;
			}
			else if (stringPosition + 1 < text.length() && text[stringPosition] == '!' && text[stringPosition + 1] == '=')
			{
				literal.is_equality = false;
				stringPosition += 2;
			}
			else
			{
				throw std::runtime_error("Syntax error: expected '=' or '!=' after attribute at position " + std::to_string(stringPosition));
			}

			if (text[stringPosition] != '"')
				throw std::runtime_error("Syntax error: expected '\"' at position " + std::to_string(stringPosition));
			stringPosition++;

			int valueStart = stringPosition;
			while (stringPosition < text.length() && text[stringPosition] != '"')
				stringPosition++;

			if (text[stringPosition] != '"')
				throw std::runtime_error("Syntax error: expected '\"' at position " + std::to_string(stringPosition));

			std::string tempValue = text.substr(valueStart, stringPosition - valueStart);

			if (corpus.string2index.find(tempValue) == corpus.string2index.end())
			{
				corpus.index2string.push_back(tempValue);
				corpus.string2index.insert(std::make_pair(tempValue, corpus.index2string.size() - 1));
			}
			literal.value = corpus.string2index.find(tempValue)->second;

			stringPosition++;
			clause.push_back(literal);
		}

		if (text[stringPosition] == ']')
			stringPosition++;
		query.push_back(clause);
	}
	return query;
}

std::vector<Match> match(const Corpus &corpus, const Query &query)
{
	std::vector<Match> matches;
	Match match;

	for (size_t sentence_index = 0; sentence_index < corpus.sentences.size(); sentence_index++)
	{
		// Get the start and end positions of the current sentence in the tokens vector
		size_t sentence_start = corpus.sentences[sentence_index];
		size_t sentence_end = corpus.sentences[sentence_index + 1];

		size_t clause_index = 0;

		for (size_t token_index = sentence_start; token_index < sentence_end; token_index++)
		{
			Clause current_clause = query[clause_index];

			for (size_t literal_index = 0; literal_index < current_clause.size(); literal_index++)
			{
				Literal current_literal = current_clause[literal_index];
				Token current_token = corpus.tokens[token_index];

				if (current_literal.attribute == "pos")
				{
					if (!(current_literal.is_equality ? current_token.pos == current_literal.value : current_token.pos != current_literal.value))
					{
						clause_index = 0;
						break;
					}
				}
				else if (current_literal.attribute == "c5")
				{
					if (!(current_literal.is_equality ? current_token.c5 == current_literal.value : current_token.c5 != current_literal.value))
					{
						clause_index = 0;

						break;
					}
				}
				else if (current_literal.attribute == "lemma")
				{
					if (!(current_literal.is_equality ? current_token.lemma == current_literal.value : current_token.lemma != current_literal.value))
					{
						clause_index = 0;
						break;
					}
				}
				else if (current_literal.attribute == "word")
				{
					if (!(current_literal.is_equality ? current_token.word == current_literal.value : current_token.word != current_literal.value))
					{
						clause_index = 0;
						break;
					}
				}
				if (literal_index + 1 == current_clause.size()) // all literal objects found in token, go to next clause
				{
					clause_index++;
				}
			}

			if (clause_index == query.size())
			{
				match.sentence = sentence_index;
				match.pos = token_index - sentence_start - (query.size() - 1);
				match.len = query.size();
				matches.push_back(match);

				clause_index = 0;
			}
		}
	}

	return matches;
}

std::vector<Match> match(Corpus &corpus, const std::string &query_string)
{
	Query query = parse_query(query_string, corpus);
	return match(corpus, query);
}

IndexSet index_lookup(const Corpus &corpus, const std::string &attribute, uint32_t value)
{
	const Index *index = nullptr;
	uint32_t Token::*attribute_p = nullptr; 

	if (attribute == "word")
	{
		index = &corpus.word_index;
		attribute_p = &Token::word;
	}
	else if (attribute == "c5")
	{
		index = &corpus.c5_index;
		attribute_p = &Token::c5;
	}
	else if (attribute == "lemma")
	{
		index = &corpus.lemma_index;
		attribute_p = &Token::lemma;
	}
	else if (attribute == "pos")
	{
		index = &corpus.pos_index;
		attribute_p = &Token::pos;
	}
	else
		throw std::runtime_error("Invalid attribute: " + attribute);

	auto lower = std::lower_bound(index->begin(), index->end(), value,
				      [&](int index_pos, uint32_t val)
				      {
					      return corpus.tokens[index_pos].*attribute_p < val;
				      });

	auto upper = std::upper_bound(lower, index->end(), value,
				      [&](uint32_t val, int index_pos)
				      {
					      return val < corpus.tokens[index_pos].*attribute_p;
				      });

	if (lower == upper)
	{
		return IndexSet{std::span<const int>(), 0}; 
	}


	std::span<const int> result_span(&(*lower), std::distance(lower, upper));
	return IndexSet{result_span, 0}; // The shift is initially set to 0
}

void build_indices(Corpus &corpus)
{
	corpus.lemma_index = build_index(corpus.tokens, &Token::lemma);
	corpus.word_index = build_index(corpus.tokens, &Token::word);
	corpus.c5_index = build_index(corpus.tokens, &Token::c5);
	corpus.pos_index = build_index(corpus.tokens, &Token::pos);
}

Index build_index(const std::vector<Token> &tokens, uint32_t Token::*attribute)
{
	Index index(tokens.size() - 1);
	for (size_t i = 0; i < tokens.size(); i++)
	{
		index[i] = i; // Initially, the index points to each token's original position
	}

	std::stable_sort(index.begin(), index.end(),
			 [&tokens, attribute](int a, int b)
			 {
				 return tokens[a].*attribute < tokens[b].*attribute;
			 });

	return index; 
}
std::vector<Match> match_single(const Corpus &corpus, const std::string &attr, const std::string &value)
{
	std::vector<Match> matches;
	Match match;
	auto it = corpus.string2index.find(value);
	if (it == corpus.string2index.end())
		return matches;

	uint32_t value_num = it->second; 

	IndexSet set = index_lookup(corpus, attr, value_num);

	if (set.elems.empty())
		return matches;

	
	for (const int token_position : set.elems)
	{
		auto sentence = std::upper_bound(corpus.sentences.begin(), corpus.sentences.end(), token_position);
		int sentence_index = std::distance(corpus.sentences.begin(), sentence) - 1;

		match.pos = token_position;	 // Token position is the starting index of the match
		match.len = 1;			 // Length is 1 since we are matching a single token
		match.sentence = sentence_index; 
		matches.push_back(match);
	}

	return matches;
}

MatchSet intersection(const MatchSet &A, const MatchSet &B)
{
	MatchSet m_set;

	if (A.complement && B.complement)
	{

		m_set.set = std::visit([](auto &&a, auto &&b) -> std::variant<DenseSet, IndexSet, ExplicitSet>
				       { return intersect(a, b); }, A.set, B.set);
		m_set.complement = true;
	}
	else if (A.complement)
	{
		m_set.set = std::visit([](auto &&a, auto &&b) -> std::variant<DenseSet, IndexSet, ExplicitSet>
				       { return difference(b, a); }, A.set, B.set);
		m_set.complement = false;
	}
	else if (B.complement)
	{

		m_set.set = std::visit([](auto &&a, auto &&b) -> std::variant<DenseSet, IndexSet, ExplicitSet>
				       { return difference(a, b); }, A.set, B.set);
		m_set.complement = false;
	}
	else
	{

		m_set.set = std::visit([](auto &&a, auto &&b) -> std::variant<DenseSet, IndexSet, ExplicitSet>
				       { return intersect(a, b); }, A.set, B.set);
		m_set.complement = false;
	}
	return m_set;
}

ExplicitSet intersect_binary(const IndexSet &A, const IndexSet &B)
{
	ExplicitSet e_set;
	for (const auto &elem : A.elems)
	{
		size_t elem_s = elem - A.shift; 

		if (std::binary_search(B.elems.begin(), B.elems.end(), elem_s + B.shift))
		{
			e_set.elems.push_back(elem_s);
		}
	}

	return e_set;
}
ExplicitSet intersect_binary(const IndexSet &A, const ExplicitSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{
		size_t elem_s = elem - A.shift;

		if (std::binary_search(B.elems.begin(), B.elems.end(), elem_s))
		{
			e_set.elems.push_back(elem_s);
		}
	}

	return e_set;
}

ExplicitSet intersect_binary(const ExplicitSet &A, const ExplicitSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{
		if (std::binary_search(B.elems.begin(), B.elems.end(), elem))
		{
			e_set.elems.push_back(elem);
		}
	}

	return e_set;
}

ExplicitSet difference_binary(const ExplicitSet &A, const ExplicitSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{

		if (!std::binary_search(B.elems.begin(), B.elems.end(), elem))
		{
			e_set.elems.push_back(elem); // Include if not found
		}
	}

	return e_set;
}
ExplicitSet difference_binary(const IndexSet &A, const IndexSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{
		size_t elem_s = elem - A.shift;

		if (!std::binary_search(B.elems.begin(), B.elems.end(), elem_s + B.shift))
		{
			e_set.elems.push_back(elem_s); // Include if not found
		}
	}

	return e_set;
}

ExplicitSet difference_binary(const IndexSet &A, const ExplicitSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{
		size_t elem_s = elem - A.shift;

		if (!std::binary_search(B.elems.begin(), B.elems.end(), elem_s))
		{
			e_set.elems.push_back(elem_s); // Include if not found
		}
	}

	return e_set;
}

ExplicitSet difference_binary(const ExplicitSet &A, const IndexSet &B)
{
	ExplicitSet e_set;

	for (const auto &elem : A.elems)
	{

		if (!std::binary_search(B.elems.begin(), B.elems.end(), elem + B.shift))
		{
			e_set.elems.push_back(elem); // Include if not found
		}
	}

	return e_set;
}

DenseSet intersect(const DenseSet &A, const DenseSet &B) 
{
	size_t first = std::max(A.first, B.first);
	size_t last = std::min(A.last, B.last);

	if (first < last)
		return DenseSet(first, last);
	else
		return DenseSet(0, 0);
}

ExplicitSet intersect(const IndexSet &A, const DenseSet &B) 
{
	ExplicitSet e_set;
	for (int elem : A.elems)
	{
		int elem_s = elem - A.shift;
		if (elem_s < B.last && elem_s >= B.first)
			e_set.elems.push_back(elem_s);
	}
	return e_set;
}
ExplicitSet intersect(const DenseSet &A, const ExplicitSet &B) // D
{
	ExplicitSet e_set;
	for (int elem : B.elems)
	{
		if (elem < A.last && elem >= A.first)
			e_set.elems.push_back(elem);
	}
	return e_set;
}
ExplicitSet intersect(const IndexSet &A, const IndexSet &B) // D
{
	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return intersect_binary(A, B);

	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		size_t p_elem = A.elems[p] - A.shift;
		size_t q_elem = B.elems[q] - B.shift;
		if (p_elem < q_elem)
			p++;
		else if (p_elem > q_elem)
			q++;
		else
		{
			e_set.elems.push_back(p_elem);
			p++;
			q++;
		}
	}
	return e_set;
}
ExplicitSet intersect(const IndexSet &A, const ExplicitSet &B) // D
{
	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return intersect_binary(A, B);

	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		int elem_s = A.elems[p] - A.shift;
		if (elem_s < B.elems[q])
			p++;

		else if (elem_s > B.elems[q])
			q++;

		else
		{
			e_set.elems.push_back(elem_s);
			p++;
			q++;
		}
	}
	
	return e_set;
}


ExplicitSet intersect(const ExplicitSet &A, const ExplicitSet &B) // D
{
	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return intersect_binary(A, B);

	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		if (A.elems[p] < B.elems[q])
			p++;
		else if (A.elems[p] > B.elems[q])
			q++;
		else
		{
			e_set.elems.push_back(A.elems[p]);
			p++;
			q++;
		}
	}
	return e_set;
}

ExplicitSet intersect(const DenseSet &A, const IndexSet &B)
{
	return intersect(B, A);
}

ExplicitSet intersect(const ExplicitSet &A, const IndexSet &B)
{
	return intersect(B, A);
}
ExplicitSet intersect(const ExplicitSet &A, const DenseSet &B)
{
	return intersect(B, A);
}

ExplicitSet difference(const ExplicitSet &A, const ExplicitSet &B)
{

	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return difference_binary(A, B);
	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		if (A.elems[p] < B.elems[q])
		{
			e_set.elems.push_back(A.elems[p]);
			p++;
		}
		else if (B.elems[q] < A.elems[p])
			q++;

		else
		{
			p++;
			q++;
		}
	}

	return e_set;
}

ExplicitSet difference(const IndexSet &A, const IndexSet &B)
{

	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return difference_binary(A, B);

	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		size_t elem_a_s = A.elems[p] - A.shift;
		size_t elem_b_s = B.elems[q] - B.shift;

		if (elem_a_s < elem_b_s)
		{
			e_set.elems.push_back(elem_a_s);
			p++;
		}
		else if (elem_a_s > elem_b_s)
			q++;
		else
		{
			p++;
			q++;
		}
	}
	while (p < A.elems.size())
	{
		e_set.elems.push_back(A.elems[p]);
		p++;
	}

	return e_set;
}

ExplicitSet difference(const IndexSet &A, const ExplicitSet &B)
{

	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return difference_binary(A, B);
	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		int elem_s = A.elems[p] - A.shift;
		if (elem_s < B.elems[q])
		{
			e_set.elems.push_back(elem_s);
			p++;
		}
		else if (B.elems[q] < elem_s)
		{
			q++;
		}
		else
		{
			p++;
			q++;
		}
	}

	return e_set;
}

ExplicitSet difference(const ExplicitSet &A, const IndexSet &B)
{

	ExplicitSet e_set;

	if (A.elems.size() * 10 < B.elems.size())
		return difference_binary(A, B);

	size_t p = 0;
	size_t q = 0;

	while (p < A.elems.size() && q < B.elems.size())
	{
		int elem_s = B.elems[q] - B.shift;

		if (A.elems[p] < elem_s)
		{
			e_set.elems.push_back(A.elems[p]);
			p++;
		}
		else if (elem_s < A.elems[p])
		{
			q++;
		}
		else
		{
			p++;
			q++;
		}
	}

	return e_set;
}

DenseSet difference(const DenseSet &A, const DenseSet &B)
{

	if (A.last <= B.first || A.first >= B.last)
	{
		return DenseSet(0, 0);
	}

	if (A.first < B.first)
	{
		return DenseSet(A.first, B.first);
	}

	return DenseSet(0, 0);
}

ExplicitSet difference(const DenseSet &A, const ExplicitSet &B)
{

	ExplicitSet e_set;

	for (int elem : B.elems)
	{
		if (elem < A.first || elem >= A.last)
		{
			e_set.elems.push_back(elem);
		}
	}

	return e_set;
}

ExplicitSet difference(const ExplicitSet &A, const DenseSet &B)
{

	ExplicitSet result;

	for (int elem : A.elems)
	{
		if (elem < B.first || elem >= B.last)
		{
			result.elems.push_back(elem);
		}
	}

	return result;
}

ExplicitSet difference(const IndexSet &A, const DenseSet &B)
{

	ExplicitSet e_set;

	for (int elem : A.elems)
	{
		int elem_s = elem - A.shift;

		if (elem_s < B.first || elem_s >= B.last)
		{
			e_set.elems.push_back(elem);
		}
	}

	return e_set;
}


ExplicitSet difference(const DenseSet &A, const IndexSet &B)
{

	ExplicitSet result;

	int p = A.first;
	size_t q = 0;

	while (p < A.last && q < B.elems.size())
	{
		int elem_s = B.elems[q] - B.shift;

		if (p < elem_s)
		{
			result.elems.push_back(p);
			p++;
		}
		else if (elem_s < p)
		{
			q++;
		}
		else
		{
			p++;
			q++;
		}
	}
	while (p <= A.last)
	{
		result.elems.push_back(p);
		p++;
	}

	return result;
}


void split_matches(std::vector<MatchSet> &match_sets, std::vector<MatchSet> &d_sets, std::vector<MatchSet> &ie_sets)
{
	for (MatchSet m_set : match_sets)
	{
		if (std::holds_alternative<DenseSet>(m_set.set))
			d_sets.push_back(m_set);

		else
			ie_sets.push_back(m_set);
	}
}

MatchSet collapse_densets(std::vector<MatchSet> &d_sets)
{
	MatchSet d_set = d_sets[0];
	for (size_t d_index = 1; d_index < d_sets.size(); d_index++)
	{
		d_set = intersection(d_set, d_sets[d_index]);
	}
	return d_set;
}
void sort_matches(std::vector<MatchSet> &matches)
{
    std::sort(matches.begin(), matches.end(),
              [](const MatchSet &a, const MatchSet &b) -> bool
              {
                  return std::visit([](const auto &setA, const auto &setB) -> bool
                                    {
                      // If both sets are DenseSet, compare by size (last - first)
                      if constexpr (std::is_same_v<std::decay_t<decltype(setA)>, DenseSet> &&
                                    std::is_same_v<std::decay_t<decltype(setB)>, DenseSet>) {
                          return (setA.last - setA.first) < (setB.last - setB.first);
                      }
                      // If both sets are ExplicitSet or IndexSet, compare by element size
                      else if constexpr ((std::is_same_v<std::decay_t<decltype(setA)>, ExplicitSet> ||
                                          std::is_same_v<std::decay_t<decltype(setA)>, IndexSet>) &&
                                         (std::is_same_v<std::decay_t<decltype(setB)>, ExplicitSet> ||
                                          std::is_same_v<std::decay_t<decltype(setB)>, IndexSet>)) {
                          return setA.elems.size() < setB.elems.size();
                      }
                      // If setA is DenseSet and setB is ExplicitSet or IndexSet, compare by size of setA
                      else if constexpr (std::is_same_v<std::decay_t<decltype(setA)>, DenseSet>) {
                          return static_cast<size_t>(setA.last - setA.first) < setB.elems.size();
                      }
                      // If setB is DenseSet and setA is ExplicitSet or IndexSet, compare by size of setB
                      else if constexpr (std::is_same_v<std::decay_t<decltype(setB)>, DenseSet>) {
                          return setA.elems.size() < static_cast<size_t>(setB.last - setB.first);
                      }

                      return false; // Default case, if no condition is met
                  }, a.set, b.set);
              });
}


MatchSet match_set(const Corpus &corpus, const Literal &literal, int shift)
{
	MatchSet m_set;
	IndexSet i_set;

	i_set = index_lookup(corpus, literal.attribute, literal.value);

	i_set.shift = shift;
	m_set.set = i_set;
	m_set.complement = !literal.is_equality;


	return m_set;
}

void match_set(const Corpus &corpus, const Clause &clause, int shift, std::vector<MatchSet> &sets)
{
	MatchSet m_set;
	if (clause.size() == 0)
	{
		m_set.set = DenseSet(0, corpus.tokens.size() - 1);
		sets.push_back(m_set);
	}

	for (size_t literal_index = 0; literal_index < clause.size(); literal_index++)
	{
		m_set = match_set(corpus, clause[literal_index], shift);
		sets.push_back(m_set);
	}

}

MatchSet match_set(const Corpus &corpus, const Query &query)
{
	MatchSet m_set;
	MatchSet d_set;
	std::vector<MatchSet> match_sets;
	std::vector<MatchSet> d_sets;
	std::vector<MatchSet> ie_sets;
	bool d_found = true;

	int shift = 0;

	for (size_t query_index = 0; query_index < query.size(); query_index++)
	{
		match_set(corpus, query[query_index], shift, match_sets);
		shift++;
	}

	split_matches(match_sets, d_sets, ie_sets);
	sort_matches(ie_sets);

	if (d_sets.size() != 0) // densesets found
		d_set = collapse_densets(d_sets);
	else
		d_found = false;


	for (size_t m_index = 0; m_index < ie_sets.size(); m_index++)
	{
		if (m_index == 0)
			m_set = ie_sets[m_index]; 
		else
			m_set = intersection(m_set, ie_sets[m_index]); 
	}
	if (d_found && ie_sets.size() != 0)
		m_set = intersection(m_set, d_set); 
	else if (ie_sets.size() == 0)
		m_set = d_set;

	if (m_set.complement)
	{
		MatchSet empty_m_set;
		empty_m_set.set = DenseSet(0, corpus.tokens.size() - 1);
		empty_m_set.complement = false;

		return intersection(empty_m_set, m_set);
	}

	return m_set;
}

std::vector<Match> match2(const Corpus &corpus, const Query &query)
{
	MatchSet m_set = match_set(corpus, query);
	std::vector<Match> matches;

	if (std::holds_alternative<DenseSet>(m_set.set))
	{
		DenseSet d_set = std::get<DenseSet>(m_set.set);

		for (int set_index = d_set.first; set_index < d_set.last; set_index++)
		{
			add_match_if_valid(matches, corpus, set_index, query.size());
		}
	}
	else if (std::holds_alternative<ExplicitSet>(m_set.set))
	{
		ExplicitSet e_set = std::get<ExplicitSet>(m_set.set);
		
		for (size_t elem_pos : e_set.elems)
		{

			add_match_if_valid(matches, corpus, elem_pos, query.size());
		}
	}
	else if (std::holds_alternative<IndexSet>(m_set.set))
	{
		IndexSet i_set = std::get<IndexSet>(m_set.set);

		for (size_t i = 0; i < i_set.elems.size(); i++)
		{
			size_t set_index = i_set.elems[i] - i_set.shift;
			add_match_if_valid(matches, corpus, set_index, query.size());
		}
	}

	return matches;
}
void add_match_if_valid(std::vector<Match> &matches, const Corpus &corpus, int set_index, int query_size)
{
	Match match;
	match.len = query_size;

	auto it = std::upper_bound(corpus.sentences.begin(), corpus.sentences.end(), set_index);
	match.pos = set_index - corpus.sentences[std::distance(corpus.sentences.begin(), it)-1];
	auto sentence = std::upper_bound(corpus.sentences.begin(), corpus.sentences.end(), set_index);
	match.sentence = std::distance(corpus.sentences.begin(), sentence) - 1;

	if (set_index + match.len <= *sentence)
	{
		matches.push_back(match);
	}
}