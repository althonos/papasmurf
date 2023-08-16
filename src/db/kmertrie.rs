

#[derive(Default, Debug, Clone)]
struct TrieNode {
    a: Option<Box<TrieNode>>,
    c: Option<Box<TrieNode>>,
    g: Option<Box<TrieNode>>,
    t: Option<Box<TrieNode>>,
}

#[derive(Default, Debug, Clone)]
pub struct KmerTrie {
    root: TrieNode,
    k: usize,
    pub nodes: usize,
}

impl KmerTrie {
    pub fn new(k: usize) -> Self {
        Self {
            k,
            root: TrieNode::default(),
            nodes: 1,
        }
    }

    pub fn insert(&mut self, word: &str) {
        assert_eq!(word.len(), self.k);
        
        let mut current_node = &mut self.root;
        for c in word.chars() {
            current_node = match c {
                'A' => current_node.a.get_or_insert_with(|| {
                    self.nodes += 1;
                    Default::default()
                }),
                'C' => current_node.c.get_or_insert_with(|| {
                    self.nodes += 1;
                    Default::default()
                }),
                'G' => current_node.g.get_or_insert_with(|| {
                    self.nodes += 1;
                    Default::default()
                }),
                'T' => current_node.t.get_or_insert_with(|| {
                    self.nodes += 1;
                    Default::default()
                }),
                _ => unreachable!(),
            };
        }
    }

    pub fn contains(&self, word: &str) -> bool {
        assert_eq!(word.len(), self.k);
        
        let mut current_node = &self.root;
        for c in word.chars() {
            let option = match c {
                'A' => &current_node.a,
                'C' => &current_node.c,
                'G' => &current_node.g,
                'T' => &current_node.t,
                _ => unreachable!(),
            };
            match option {
                Some(node) => current_node = node,
                None => return false,
            }
        }

        true
    }

    pub fn fuzzy_search(&self, word: &str, max_mismatches: usize) -> Vec<(String, usize)> {
        // implement a depth-first search algorithm

        struct State<'a> {
            node: &'a TrieNode,
            mismatches: usize,
            position: usize,
            prefix: String,
        }

        let mut results = Vec::new();

        let mut bytes = word.as_bytes();
        let mut stack = Vec::new();
        stack.push(State{
            node: &self.root, 
            mismatches: 0, 
            position: 0,
            prefix: String::new(),
        });

        while let Some(state) = stack.pop() {
            if state.position == self.k {
                results.push((state.prefix, state.mismatches));
            } else {
                let c = bytes[state.position];
                if let Some(node) = &state.node.a {
                    if c == b'A' || state.mismatches < max_mismatches {
                        stack.push(State {
                            node,
                            mismatches: state.mismatches + (c != b'A') as usize,
                            position: state.position + 1,
                            prefix: format!("{}A", &state.prefix),
                        })
                    }
                }
                if let Some(node) = &state.node.c {
                    if c == b'C' || state.mismatches < max_mismatches {
                        stack.push(State {
                            node,
                            mismatches: state.mismatches + (c != b'C') as usize,
                            position: state.position + 1,
                            prefix: format!("{}C", &state.prefix),
                        })
                    }
                }
                if let Some(node) = &state.node.g {
                    if c == b'G' || state.mismatches < max_mismatches {
                        stack.push(State {
                            node,
                            mismatches: state.mismatches + (c != b'G') as usize,
                            position: state.position + 1,
                            prefix: format!("{}G", &state.prefix),
                        })
                    }
                }
                if let Some(node) = &state.node.t {
                    if c == b'T' || state.mismatches < max_mismatches {
                        stack.push(State {
                            node,
                            mismatches: state.mismatches + (c != b'T') as usize,
                            position: state.position + 1,
                            prefix: format!("{}T", &state.prefix),
                        })
                    }
                }
            }
        }

        results
        // let mut  = &self.root;

        // for c in word.chars() {

        //     let option = match c {
        //         'A' => &current_node.a,
        //         'C' => &current_node.c,
        //         'G' => &current_node.g,
        //         'T' => &current_node.t,
        //         _ => unreachable!(),
        //     };

        //     match option {
        //         Some(node) => current_node = node,
        //         None => return false,
        //     }
        // }

        // current_node.is_end_of_word
        


    }
}

#[cfg(test)]
mod test {

    use super::*;

    #[test]
    fn test_contains() {
        let mut trie = KmerTrie::new(4);
        assert!(!trie.contains("ATTA"));
        assert!(!trie.contains("ATTT"));

        trie.insert("ATTA");
        trie.insert("ATGC");
        assert!(!trie.contains("ATTC"));
        assert!(!trie.contains("CGTC"));
        assert!(trie.contains("ATTA"));
        assert!(trie.contains("ATGC"));
    }

    #[test]
    fn test_fuzzy_search() {
        let mut trie = KmerTrie::new(4);
        trie.insert("ATTA");
        trie.insert("ATTT");
        trie.insert("ATGC");

        let mut results = trie.fuzzy_search("ATTA", 1);
        println!("{:?}", results);
        results.sort();

        assert_eq!(results.len(), 2);
        assert_eq!(results[0], "ATTA");
        assert_eq!(results[1], "ATTT");
    }
}