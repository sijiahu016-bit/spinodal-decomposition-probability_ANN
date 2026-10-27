function z = predict(net,x)
probabilities = net(x');
a = max(probabilities,[],1);
z = a(:);
end