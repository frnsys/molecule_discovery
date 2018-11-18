import React, {Component} from 'react';
import {render} from 'react-dom';

function get(path, onSuccess) {
  fetch(path, {
    headers: {
      'Accept': 'application/json',
      'Content-Type': 'application/json'
    },
    method: 'get'
  })
    .then(res => res.json())
    .then((data) => onSuccess(data))
    .catch(err => { throw err });
}

const Compound = (props) => {
  var formula = props.formula.match(/(\d+|[^\d]+)/g).map((f, i) => isNaN(f) ? f : <sub key={i}>{f}</sub>);
  return (
    <li className='compound'>
      <figure>
        <img src={`img/${props.id}.png`} />
      </figure>
      <div className='meta'>
        <div>{formula}</div>
        <div><span>ID:</span> {props.id}</div>
        <div><span>SMILES:</span> {props.smiles}</div>
        <div><span>BIRTHDAY:</span> {props.created_at}</div>
        <div><span>ATC:</span> {props.atc_code}
          <ul className='atc_code'>
            {atc_descs[props.atc_code].slice(0, 3).map((desc, i) => {
              return <li key={i}>{desc}</li>;
            })}
          </ul>
        </div>
      </div>
    </li>
  );
}

class Cluster extends Component {
  constructor(props) {
    super(props);
    this.state = {
      open: false
    };
  }

  toggleDetails() {
    if (!this.state.open) {
      if (this.state.members) {
        this.setState({open: true});
      } else {
        get(`clusters/${this.props.label}.json`, (data) => {
          this.setState({open: true, ...data});
        });
      }
    } else {
      this.setState({open: false});
    }
  }

  render() {
    return (
      <div open={this.state.open}>
        <header onClick={this.toggleDetails.bind(this)}>
          <h2>{this.props.label.replace('__', ': ')} <span>{this.props.size}</span></h2>
          <h3>{this.props.name}</h3>
        </header>
        {this.state.open && this.state.members &&
          <div>
            <ul className='atc_codes group'>
              {Object.keys(this.state.atc_codes).sort().map((code) => {
                return <li key={code}>
                  {`${code} ${this.state.atc_codes[code]}`}
                  <ul className='meta tip'>
                    {atc_descs[code].map((desc, i) => {
                      return <li key={i}>{desc}</li>;
                    })}
                  </ul>
                </li>;
              })}
            </ul>
            <ul className='group'>
              {this.state.members.map((c, i) => {
                return <Compound key={i} {...c} />;
              })}
            </ul>
          </div>}
      </div>
    );
  }
}

class App extends Component {
  render() {
    return (
      <div>
        {this.props.clusters.map((clus) => {
          return <Cluster key={clus.label} {...clus} />;
        })}
      </div>
    );
  }
}

// Sticky cluster headers
let stickyHeader = null;
let stickyHeaderEl = document.getElementById('sticky');
window.onscroll = () => {
  let headers = Array.from(document.getElementsByTagName('header'));
  let closest = headers.filter((el) => el.getBoundingClientRect().top < 0).slice(-1)[0];
  if (closest && closest.parentElement.getAttribute('open') !== null) {
    if (stickyHeader !== closest) {
      stickyHeader = closest;
      while (stickyHeaderEl.firstChild) {
        stickyHeaderEl.removeChild(stickyHeaderEl.firstChild);
      }
      stickyHeaderEl.appendChild(stickyHeader.cloneNode(true));
    }
    stickyHeaderEl.style.display = 'block';
  } else {
    stickyHeaderEl.style.display = 'none';
  }
}

let atc_descs = {};
get(`data/atc_descs.json`, (descs) => {
  atc_descs = descs;
  get(`data/clusters.json`, (clusters) => {
    let main = document.getElementById('main');
    render(<App clusters={clusters}/>, main);
  });
});
